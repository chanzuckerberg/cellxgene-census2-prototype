#!/usr/bin/env node
// Refresh the tracker fields of snapshot.json from the Linear API.
// Needs LINEAR_API_KEY: a personal API key with read access to the AI Products team.
// Keeps the dashboard's own fields (iteration texts and statuses, workstream names, colours).
// The workflow runs it before the build on the schedule; the committed snapshot is the baseline.
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const here = dirname(fileURLToPath(import.meta.url));
const snapshotPath = join(here, 'snapshot.json');
const PROJECT = '95fcb2ae-ebee-4e80-8bae-a348e4e1899d'; // Census on Zarr, P-AIP-165
const QUERY = `query Project($id: String!, $after: String) {
  project(id: $id) {
    projectMilestones { nodes { name targetDate } }
    issues(first: 100, after: $after) {
      pageInfo { hasNextPage endCursor }
      nodes { identifier title url dueDate description state { name } assignee { name } parent { identifier } projectMilestone { name } }
    }
  }
}`;

export async function fetchProject(token, fetchImpl = fetch) {
  const issues = [];
  let milestones = [];
  let after = null;
  do {
    let response;
    for (let attempt = 0; ; attempt++) {
      response = await fetchImpl('https://api.linear.app/graphql', {
        method: 'POST',
        headers: { Authorization: token, 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: QUERY, variables: { id: PROJECT, after } }),
      });
      if (response.ok) break;
      const body = await response.text();
      // Linear reports a rate limit as HTTP 400 with the error code RATELIMITED.
      const limited = response.status === 429 || (response.status === 400 && body.includes('RATELIMITED'));
      if (!(limited || response.status >= 500) || attempt >= 5) throw new Error(`Linear ${response.status}: ${body}`);
      const retryAfter = response.headers.get('retry-after');
      const reset = Number(response.headers.get('x-ratelimit-requests-reset'));
      const wait = retryAfter !== null ? 1000 * Number(retryAfter)
        : reset > Date.now() ? Math.min(reset - Date.now(), 60_000) : 1000 * 2 ** attempt;
      await new Promise((resolve) => setTimeout(resolve, wait));
    }
    const { data, errors } = await response.json();
    if (errors?.length) throw new Error(`Linear: ${errors.map((error) => error.message).join('; ')}`);
    if (!data?.project) throw new Error('Linear: the key cannot read the Census on Zarr project.');
    milestones = data.project.projectMilestones.nodes;
    issues.push(...data.project.issues.nodes);
    const { hasNextPage, endCursor } = data.project.issues.pageInfo;
    after = hasNextPage ? endCursor : null;
  } while (after);
  return { milestones, issues };
}

const numericId = (title) => title.match(/^(?:Phase\s+)?(\d+(?:\.\d+)?)/)?.[1];
const byId = (a, b) => a.id.split('.').map(Number).reduce((cmp, n, i) => cmp || n - Number(b.id.split('.')[i] ?? -1), 0);
const iterationId = (name) => name?.match(/^Iteration (\d+):/)?.[1];
const plannedWindow = (description) => description?.match(/\*\*Planned window:\*\*\s*(\d{4}-\d{2}-\d{2})\s*→\s*(\d{4}-\d{2}-\d{2})/);
// Mike is not in Linear yet, so his tasks name him in an "Owner:" or "Assignee:" line instead.
const owner = (issue) => issue.assignee?.name
  ?? issue.description?.match(/^\*\*(?:Owner|Assignee):\*\*[ \t]*([^(\r\n]*[^(\s])/m)?.[1]
  ?? 'Unassigned';

export function mergeSnapshot(snapshot, { milestones, issues }) {
  const order = snapshot.iterations.map((item) => item.id);
  const dates = new Map(milestones.filter((m) => iterationId(m.name)).map((m) => [iterationId(m.name), m.targetDate ?? null]));
  const iterations = snapshot.iterations.map((item) => (dates.has(item.id) ? { ...item, date: dates.get(item.id) } : item));
  const keep = Object.fromEntries(snapshot.epics.map((epic) => [epic.id, epic]));
  const epics = issues.filter((issue) => !issue.parent && /^Phase \d+:/.test(issue.title)).map((issue) => {
    const id = numericId(issue.title);
    const prior = keep[id] ?? {};
    const window = plannedWindow(issue.description);
    return {
      id, title: issue.title,
      color: prior.color ?? '#607d8b',
      url: issue.url,
      status: issue.state.name,
      iterations: [...(prior.iterations ?? [])],
      trackerStart: window?.[1] ?? null,
      trackerEnd: window?.[2] ?? issue.dueDate ?? null,
      workstream: prior.workstream ?? issue.title.replace(/^Phase \d+:\s*/, ''),
      identifier: issue.identifier,
    };
  }).sort(byId);
  const epicByIdentifier = Object.fromEntries(epics.map((epic) => [epic.identifier, epic]));
  const tasks = issues.map((issue) => {
    const epic = epicByIdentifier[issue.parent?.identifier];
    if (!epic) return null;
    const id = numericId(issue.title) ?? issue.title;
    return {
      id, title: issue.title,
      url: issue.url,
      status: issue.state.name,
      epic: epic.id,
      iteration: iterationId(issue.projectMilestone?.name) ?? epic.iterations[0] ?? 'later',
      owner: owner(issue),
      trackerStart: null,
      trackerEnd: issue.dueDate ?? null,
    };
  }).filter(Boolean).sort(byId);
  // The build requires each task's iteration among its epic's iterations.
  for (const task of tasks) {
    const epic = epics.find((item) => item.id === task.epic);
    if (!epic.iterations.includes(task.iteration)) epic.iterations.push(task.iteration);
  }
  epics.forEach((epic) => {
    epic.iterations.sort((a, b) => order.indexOf(a) - order.indexOf(b));
    delete epic.identifier;
  });
  return { ...snapshot, capturedAt: new Date().toISOString().replace(/\.\d{3}Z$/, 'Z'), iterations, epics, tasks };
}

export async function refresh(token, fetchImpl = fetch) {
  const snapshot = JSON.parse(readFileSync(snapshotPath, 'utf8'));
  const next = mergeSnapshot(snapshot, await fetchProject(token, fetchImpl));
  writeFileSync(snapshotPath, `${JSON.stringify(next, null, 2)}\n`);
  return next;
}

if (process.argv[1] && fileURLToPath(import.meta.url) === process.argv[1]) {
  const token = process.env.LINEAR_API_KEY;
  if (!token) { console.error('LINEAR_API_KEY is not set.'); process.exit(1); }
  const next = await refresh(token);
  console.log(`Refreshed snapshot.json: ${next.epics.length} epics, ${next.tasks.length} tasks.`);
}
