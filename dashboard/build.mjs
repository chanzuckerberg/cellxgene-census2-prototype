import assert from 'node:assert/strict';
import { readFile, writeFile } from 'node:fs/promises';

const read = (name) => readFile(new URL(name, import.meta.url), 'utf8');
const [template, css, script, raw, paletteRaw, favicon] = await Promise.all([
  read('template.html'), read('style.css'), read('app.js'), read('snapshot.json'), read('palettes.json'), read('favicon.svg'),
]);
const data = JSON.parse(raw);
const palettes = JSON.parse(paletteRaw);
assert.equal(data.schemaVersion, 1);
assert.ok(Number.isFinite(Date.parse(data.capturedAt)));
const unique = (items) => {
  const ids = new Set(items.map((item) => item.id));
  assert.equal(ids.size, items.length, 'IDs must be unique');
  return ids;
};
const iterations = unique(data.iterations);
const epics = unique(data.epics);
assert.equal(new Set(data.epics.map((epic) => epic.workstream)).size, data.epics.length, 'Workstream names must be unique');
for (const palette of Object.values(palettes)) {
  assert.match(palette.ink, /^#[0-9a-f]{6}$/i);
  assert.deepEqual(new Set(Object.keys(palette.colors)), epics);
  assert.equal(new Set(Object.values(palette.colors)).size, epics.size);
  Object.values(palette.colors).forEach((color) => assert.match(color, /^#[0-9a-f]{6}$/i));
}
assert.equal(new Set(data.epics.map((epic) => epic.color.toLowerCase())).size, data.epics.length, 'Epic colours must be unique');
unique(data.tasks);
const safeLink = (value) => {
  assert.match(value, /^https:\/\/app\.notion\.com\/p\/[0-9a-f]{32}$/i, 'Use the Notion page URL: https://app.notion.com/p/<page-id>');
};
for (const key of ['design', 'outline', 'project', 'plan', 'tracker', 'iterationProposal']) {
  safeLink(data.sources[key]);
}
assert.match(data.sources.twoPager, /^https:\/\/docs\.google\.com\/document\/d\/[\w-]+\/edit(?:\?usp=sharing)?$/, 'Use a Google Docs document URL');
assert.match(data.sources.slides, /^https:\/\/docs\.google\.com\/presentation\/d\/[\w-]+\/edit(?:\?usp=sharing)?$/, 'Use a Google Slides presentation URL');
for (const iteration of data.iterations) {
  assert.ok(iteration.date === null || /^\d{4}-\d{2}-\d{2}$/.test(iteration.date));
}
for (const epic of data.epics) {
  assert.ok(typeof epic.workstream === 'string' && epic.workstream.trim(), 'Each epic needs a workstream name');
  assert.match(epic.color, /^#[0-9a-f]{6}$/i);
  assert.ok(typeof epic.status === 'string' && epic.status.trim(), 'Each epic needs a tracker status');
  safeLink(epic.url);
  assert.ok(epic.iterations.every((id) => iterations.has(id)));
}
for (const task of data.tasks) {
  safeLink(task.url);
  assert.ok(epics.has(task.epic));
  assert.ok(iterations.has(task.iteration));
  assert.ok(data.epics.find((epic) => epic.id === task.epic).iterations.includes(task.iteration));
}
const serialized = JSON.stringify(data).replaceAll('<', '\\u003c');
const replacements = { CSS: css, DATA: serialized, PALETTES: JSON.stringify(palettes).replaceAll('<', '\\u003c'), SCRIPT: script, FAVICON: `data:image/svg+xml,${encodeURIComponent(favicon.trim())}` };
const html = template.replace(/\/\*__(CSS|DATA|PALETTES|SCRIPT|FAVICON)__\*\//g, (_, key) => replacements[key]);
await writeFile(new URL('index.html', import.meta.url), html);
console.log(`Built dashboard/index.html: ${data.epics.length} epics, ${data.tasks.length} tasks.`);
