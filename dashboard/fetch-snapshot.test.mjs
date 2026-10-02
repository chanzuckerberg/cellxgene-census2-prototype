import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { fetchProject, mergeSnapshot } from './fetch-snapshot.mjs';

const issue = (identifier, title, status, parent, extra = {}) => ({
  identifier, title,
  url: `https://linear.app/biohub/issue/${identifier}/slug`,
  dueDate: extra.dueDate ?? null,
  description: extra.description ?? '',
  state: { name: status },
  assignee: extra.assignee ? { name: extra.assignee } : null,
  parent: parent ? { identifier: parent } : null,
  projectMilestone: extra.milestone ? { name: extra.milestone } : null,
});

test('Linear issues replace the tracker fields and keep the dashboard fields', async () => {
  const snapshot = JSON.parse(await readFile(new URL('snapshot.json', import.meta.url), 'utf8'));
  const prior = snapshot.epics.find((epic) => epic.id === '5');
  const milestones = [
    { name: 'Iteration 1: End-to-end skeleton', targetDate: '2026-10-01' },
    { name: 'Iteration 4: Continuous ingestion and private data', targetDate: null },
    { name: 'Some other milestone', targetDate: '2026-12-31' },
  ];
  const issues = [
    issue('AIP-129', 'Phase 5: The catalog', 'Done', null, { dueDate: '2026-10-30', description: '**Planned window:** 2026-10-19 → 2026-10-30\n' }),
    issue('AIP-153', '5.3 Build the exclusion path in the SDK', 'In Progress', 'AIP-129', { assignee: 'Ronen Kalo', milestone: 'Iteration 3: Scale and layout' }),
    issue('AIP-901', '5.9 A new task', 'Backlog', 'AIP-129'),
    issue('AIP-156', '5.6 Build catalog rows from store summary records', 'Backlog', 'AIP-129', { milestone: 'Iteration 0: Design and measurement', description: '**Assignee:** Mike Lin (not yet in Linear)\n' }),
    issue('AIP-155', '5.5 Implement the publication protocol', 'Backlog', 'AIP-129', { description: '**Owner:**\nThe next paragraph.\n' }),
    issue('AIP-999', 'Orphan task', 'Backlog', 'AIP-000'),
    issue('AIP-998', 'Not a phase', 'Backlog', null),
    issue('AIP-997', 'Phase 5 follow-up', 'Backlog', null),
  ];
  const next = mergeSnapshot(snapshot, { milestones, issues });
  assert.deepEqual(next.epics, [{
    ...prior, status: 'Done', url: 'https://linear.app/biohub/issue/AIP-129/slug',
    trackerStart: '2026-10-19', trackerEnd: '2026-10-30',
    iterations: ['0', ...prior.iterations.filter((id) => id !== '0')],
  }]);
  assert.deepEqual(next.tasks.map((task) => [task.id, task.status, task.iteration, task.owner]), [
    ['5.3', 'In Progress', '3', 'Ronen Kalo'],
    // No milestone: the epic's first iteration, not the committed one.
    ['5.5', 'Backlog', prior.iterations[0], 'Unassigned'],
    ['5.6', 'Backlog', '0', 'Mike Lin'],
    ['5.9', 'Backlog', prior.iterations[0], 'Unassigned'],
  ]);
  assert.equal(next.iterations.find((item) => item.id === '1').date, '2026-10-01');
  assert.equal(next.iterations.find((item) => item.id === '4').date, null);
  assert.deepEqual(next.iterations.find((item) => item.id === '2'), snapshot.iterations.find((item) => item.id === '2'));
  assert.match(next.capturedAt, /^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$/);
});

test('fetchProject follows the issue cursor and stops on API errors', async () => {
  const pages = [
    { project: { projectMilestones: { nodes: [] }, issues: { pageInfo: { hasNextPage: true, endCursor: 'c1' }, nodes: [{ identifier: 'AIP-1' }] } } },
    { project: { projectMilestones: { nodes: [] }, issues: { pageInfo: { hasNextPage: false, endCursor: null }, nodes: [{ identifier: 'AIP-2' }] } } },
  ];
  const cursors = [];
  const fetchImpl = async (_url, init) => {
    cursors.push(JSON.parse(init.body).variables.after);
    return { ok: true, json: async () => ({ data: pages[cursors.length - 1] }) };
  };
  const result = await fetchProject('key', fetchImpl);
  assert.deepEqual(result.issues.map((item) => item.identifier), ['AIP-1', 'AIP-2']);
  assert.deepEqual(cursors, [null, 'c1']);
  const failing = async () => ({ ok: true, json: async () => ({ errors: [{ message: 'Authentication required' }] }) });
  await assert.rejects(fetchProject('key', failing), /Authentication required/);
  const hidden = async () => ({ ok: true, json: async () => ({ data: { project: null } }) });
  await assert.rejects(fetchProject('key', hidden), /cannot read/);
});

test('fetchProject retries a Linear rate limit reported as HTTP 400', async () => {
  const page = { project: { projectMilestones: { nodes: [] }, issues: { pageInfo: { hasNextPage: false, endCursor: null }, nodes: [{ identifier: 'AIP-1' }] } } };
  const headers = (values) => ({ get: (name) => values[name] ?? null });
  let calls = 0;
  const fetchImpl = async () => (++calls === 1
    ? { ok: false, status: 400, headers: headers({ 'retry-after': '0' }), text: async () => '{"errors":[{"extensions":{"code":"RATELIMITED"}}]}' }
    : { ok: true, json: async () => ({ data: page }) });
  const result = await fetchProject('key', fetchImpl);
  assert.equal(calls, 2);
  assert.deepEqual(result.issues.map((item) => item.identifier), ['AIP-1']);
  const badRequest = async () => ({ ok: false, status: 400, headers: headers({}), text: async () => 'Bad request' });
  await assert.rejects(fetchProject('key', badRequest), /Linear 400: Bad request/);
});
