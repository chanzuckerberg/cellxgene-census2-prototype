import assert from 'node:assert/strict';
import { readFile, writeFile } from 'node:fs/promises';

const read = (name) => readFile(new URL(name, import.meta.url), 'utf8');
const [template, css, script, raw] = await Promise.all([
  read('template.html'), read('style.css'), read('app.js'), read('snapshot.json'),
]);
const data = JSON.parse(raw);
assert.equal(data.schemaVersion, 1);
assert.ok(Number.isFinite(Date.parse(data.capturedAt)));
const unique = (items) => {
  const ids = new Set(items.map((item) => item.id));
  assert.equal(ids.size, items.length, 'IDs must be unique');
  return ids;
};
const iterations = unique(data.iterations);
const epics = unique(data.epics);
assert.equal(new Set(data.epics.map((epic) => epic.color.toLowerCase())).size, data.epics.length, 'Epic colours must be unique');
unique(data.tasks);
const safeLink = (value) => {
  const url = new URL(value);
  assert.equal(url.protocol, 'https:');
  assert.ok(['app.notion.com', 'www.notion.so'].includes(url.hostname));
};
Object.values(data.sources).forEach(safeLink);
for (const iteration of data.iterations) {
  assert.ok(iteration.date === null || /^\d{4}-\d{2}-\d{2}$/.test(iteration.date));
}
for (const epic of data.epics) {
  assert.match(epic.color, /^#[0-9a-f]{6}$/i);
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
const replacements = { CSS: css, DATA: serialized, SCRIPT: script };
const html = template.replace(/\/\*__(CSS|DATA|SCRIPT)__\*\//g, (_, key) => replacements[key]);
await writeFile(new URL('index.html', import.meta.url), html);
console.log(`Built dashboard/index.html: ${data.epics.length} epics, ${data.tasks.length} tasks.`);
