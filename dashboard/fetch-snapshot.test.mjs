import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { mergeSnapshot } from './fetch-snapshot.mjs';

const page = (id, title, status, parent, extra = {}) => ({
  id,
  properties: {
    Task: { title: [{ plain_text: title }] },
    Status: { status: { name: status } },
    'Parent Task': { relation: parent ? [{ id: parent }] : [] },
    'Due Date': { date: extra.date ?? null },
    Owner: { people: extra.owner ? [{ name: extra.owner }] : extra.nameless ? [{ id: 'u1' }] : [] },
  },
});

test('tracker pages replace statuses and keep the dashboard mapping', async () => {
  const snapshot = JSON.parse(await readFile(new URL('snapshot.json', import.meta.url), 'utf8'));
  const epic5 = page('e5000000-0000-0000-0000-000000000005', 'Phase 5: The catalog', 'Done', null, { date: { start: '2026-10-19', end: '2026-10-30' } });
  const tasks = [
    page('t5300000-0000-0000-0000-000000000053', '5.3 Build the exclusion path in the SDK', 'Done', epic5.id, { owner: 'Ronen Kalo' }),
    page('t5900000-0000-0000-0000-000000000059', '5.9 A new task', 'Backlog', epic5.id),
    page('t5600000-0000-0000-0000-000000000056', '5.6 Build catalog rows from store summary records', 'Backlog', epic5.id, { nameless: true }),
    page('t9999999-0000-0000-0000-000000000099', 'Orphan task', 'Backlog', 'not-an-epic'),
  ];
  const next = mergeSnapshot(snapshot, [epic5], tasks);
  const prior = snapshot.epics.find((epic) => epic.id === '5');
  assert.deepEqual(next.epics, [{ ...prior, status: 'Done', url: 'https://app.notion.com/p/e5000000000000000000000000000005' }]);
  assert.deepEqual(next.tasks.map((task) => [task.id, task.status, task.iteration, task.owner]), [
    ['5.3', 'Done', '3', 'Ronen Kalo'],
    ['5.6', 'Backlog', '3', snapshot.tasks.find((task) => task.id === '5.6').owner],
    ['5.9', 'Backlog', prior.iterations[0], 'Unassigned'],
  ]);
  assert.deepEqual(next.iterations, snapshot.iterations);
  assert.match(next.capturedAt, /^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$/);
});
