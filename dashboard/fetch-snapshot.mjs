#!/usr/bin/env node
// Refresh the tracker fields of snapshot.json from the Notion API.
// Needs NOTION_TOKEN from an internal integration connected to the tracker database.
// Keeps the dashboard's own fields (iterations, workstream names, colours, the iteration mapping).
// The workflow runs it before the build on the schedule; the committed snapshot is the baseline.
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const snapshotPath = join(here, 'snapshot.json');
const DATA_SOURCE = 'b0816a75-3789-82f1-91e3-87989b8c7954';
const PROJECT = '38f16a75-3789-806b-8fb0-e21950fd3659';

export async function queryAll(filter, token, fetchImpl = fetch) {
  const results = [];
  let cursor;
  do {
    let response;
    for (let attempt = 0; ; attempt++) {
      response = await fetchImpl(`https://api.notion.com/v1/data_sources/${DATA_SOURCE}/query`, {
        method: 'POST',
        headers: { Authorization: `Bearer ${token}`, 'Notion-Version': '2025-09-03', 'Content-Type': 'application/json' },
        body: JSON.stringify({ filter, page_size: 100, ...(cursor ? { start_cursor: cursor } : {}) }),
      });
      if (response.ok) break;
      if (![429, 529].includes(response.status) || attempt >= 5) throw new Error(`Notion ${response.status}: ${await response.text()}`);
      await new Promise((resolve) => setTimeout(resolve, 1000 * Number(response.headers.get('retry-after') ?? 2 ** attempt)));
    }
    const page = await response.json();
    results.push(...page.results);
    cursor = page.has_more ? page.next_cursor : undefined;
  } while (cursor);
  return results;
}

const text = (property) => (property?.title ?? property?.rich_text ?? []).map((part) => part.plain_text).join('');
const person = (page, fallback) => {
  // A people value can be a partial user object without a name; keep the prior owner only then.
  const people = page.properties.Owner?.people?.length ? page.properties.Owner.people : page.properties.Assignee?.people ?? [];
  return people.length ? people[0].name ?? fallback ?? 'Unassigned' : 'Unassigned';
};
const pageUrl = (page) => `https://app.notion.com/p/${page.id.replaceAll('-', '')}`;
const numericId = (title) => title.match(/^(?:Phase\s+)?(\d+(?:\.\d+)?)/)?.[1];
const byId = (a, b) => a.id.split('.').map(Number).reduce((cmp, n, i) => cmp || n - Number(b.id.split('.')[i] ?? -1), 0);

export function mergeSnapshot(snapshot, epicPages, taskPages) {
  const keep = Object.fromEntries(snapshot.epics.map((epic) => [epic.id, epic]));
  const epics = epicPages.map((page) => {
    const title = text(page.properties.Task);
    const id = numericId(title);
    const prior = keep[id] ?? {};
    return {
      id, title,
      color: prior.color ?? '#607d8b',
      url: pageUrl(page),
      status: page.properties.Status?.status?.name ?? 'Backlog',
      iterations: prior.iterations ?? [],
      trackerStart: page.properties['Due Date']?.date?.start ?? null,
      trackerEnd: page.properties['Due Date']?.date?.end ?? null,
      workstream: prior.workstream ?? title.replace(/^Phase \d+:\s*/, ''),
      pageId: page.id,
    };
  }).sort(byId);
  const epicByPageId = Object.fromEntries(epics.map((epic) => [epic.pageId, epic]));
  const keepTask = Object.fromEntries(snapshot.tasks.map((task) => [task.id, task]));
  const tasks = taskPages.map((page) => {
    const title = text(page.properties.Task);
    const epic = epicByPageId[page.properties['Parent Task']?.relation?.[0]?.id];
    if (!epic) return null;
    const id = numericId(title) ?? title;
    return {
      id, title,
      url: pageUrl(page),
      status: page.properties.Status?.status?.name ?? 'Backlog',
      epic: epic.id,
      iteration: keepTask[id]?.iteration ?? epic.iterations[0] ?? 'later',
      owner: person(page, keepTask[id]?.owner),
      trackerStart: page.properties['Due Date']?.date?.start ?? null,
      trackerEnd: page.properties['Due Date']?.date?.end ?? null,
    };
  }).filter(Boolean).sort(byId);
  epics.forEach((epic) => delete epic.pageId);
  return { ...snapshot, capturedAt: new Date().toISOString().replace(/\.\d{3}Z$/, 'Z'), epics, tasks };
}

export async function refresh(token, fetchImpl = fetch) {
  const snapshot = JSON.parse(readFileSync(snapshotPath, 'utf8'));
  const epicPages = await queryAll({ property: 'Project', relation: { contains: PROJECT } }, token, fetchImpl);
  const taskPages = await queryAll(
    { or: epicPages.map((page) => ({ property: 'Parent Task', relation: { contains: page.id } })) },
    token, fetchImpl,
  );
  const next = mergeSnapshot(snapshot, epicPages, taskPages);
  writeFileSync(snapshotPath, `${JSON.stringify(next, null, 2)}\n`);
  return next;
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const token = process.env.NOTION_TOKEN;
  if (!token) { console.error('NOTION_TOKEN is not set.'); process.exit(1); }
  const next = await refresh(token);
  console.log(`Refreshed snapshot.json: ${next.epics.length} epics, ${next.tasks.length} tasks.`);
}
