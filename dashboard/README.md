# Census 2 project dashboard

The dashboard shows the iteration gates from the presentation outline. Each iteration links to
the relevant Notion epics and tasks. The dashboard uses a dated snapshot of the tracker.

[GitHub Pages hosts the team dashboard](https://scaling-dollop-pzyy6pk.pages.github.io/)
with private access. The
[Notion dashboard page](https://app.notion.com/p/3ea16a75378981019542e9e2a3db2150)
also holds the original HTML snapshot.

Open `index.html` in a browser. The file also works as an HTML attachment in Notion.
It contains its data, styles, and script. It makes no API requests and needs no credentials.
The build embeds `favicon.svg` as the browser tab icon. The header links to the current design
record, the query diagram, Meg's 2-pager, the slides, the presentation outline, the project, and
the tracker.
Click the Census 2 header to return home. This closes the drawer, selects the iteration in progress, clears
task filters, closes the colour settings, and returns the page and tables to their start positions.
The selected theme and workstream colours stay unchanged.

Click a workstream band to highlight that workstream across the timeline. A right-side drawer shows its
saved status, related iterations, linked tasks, and an explicit link to the epic in Notion.
The timeline stays interactive on desktop. On smaller screens, the drawer opens as a modal.
Close it with the close button or Escape. Focus returns to the workstream band.
The task table also shows the workstream's linked tasks. Its search and status filters still apply.
Click an iteration label or the Back button to close the drawer and show the iteration.
A URL with `#workstream-4` opens Writer. Workstream selection does not open Notion automatically.
Each timeline row shows what the iteration delivers, in the outline's words. Select the
iteration to read the deliverable list, the gate, and the proposed demo or release summary.
These texts come from the outline and are for team review. They do not establish completion
or change the proposed gate dates.
Workstreams span iterations. Each maps to one tracker epic. The `workstream` field in the snapshot
holds its short display name, such as SDK or Catalog. The drawer retains the original epic title
and phase number for traceability. Dashboard labels do not rename records in Notion.

## Build

Run these commands from the repository root. They need Node.js 18 or later.

```sh
node dashboard/build.mjs
node --test dashboard/build.test.mjs
```

The build validates the snapshot and writes `dashboard/index.html`. Commit the generated file
with the source files so the preview opens without a build step.
Page links must use `https://app.notion.com/p/<page-id>`. The build rejects URLs without `/p/`.

## Sources and dates

- The [presentation outline](https://app.notion.com/p/3e616a753789813491d6ffa66a4c3998)
  supplies the iteration names, proposed gate dates, and outline statuses.
- The [project](https://app.notion.com/p/38f16a753789806b8fb0e21950fd3659)
  supplies the epic membership.
- The [tracker](https://app.notion.com/p/95716a7537898228b0150198eacd21dc)
  supplies task titles, statuses, owners, and parent relations.
- The [plan](https://app.notion.com/p/3b216a75378981c9b7f3d954789aefe4)
  and [iteration proposal](https://app.notion.com/p/68916a75378982a3908e018e928d7eb7)
  supply the workstream deliverables used for the proposed iteration mapping.

The tracker has no iteration property. `snapshot.json` holds the proposed mapping.
An epic can support several iterations. Each child task has one proposed iteration.
The dashboard labels this mapping for review. It does not change the tracker.

Each epic has an original colour in `snapshot.json`. Workstream names identify the chart bands.
Each band shows the epic's saved tracker status. The status applies to the whole epic, not its
completion within each iteration. The bands do not assign separate epic dates.
Workstream bars use status icons: a circled check for Done, a half-filled circle for In Progress,
an empty circle for Backlog, and an eye for In Review. The legend identifies the displayed icons.
Tooltips, accessible labels, the drawer, and task tables retain the full status text.

Midnight is the default theme. A saved choice or URL theme takes precedence.
The theme and workstream colour pickers start closed. Click the colour-settings icon to open them.
The dashboard theme selector offers Slate, Sand, and Midnight. System uses the original green theme
and follows the device's light or dark setting. The browser saves the choice when storage is
available. A URL with `?theme=slate`, `?theme=sand`, or `?theme=midnight` opens that theme.

The workstream colour selector compares the current palette with muted and vivid variants from `palettes.json`.
The selection applies to the chart, legend, and links. The browser saves the choice when storage
is available. A URL with `?palette=muted` or `?palette=vivid` opens the selected variant.
Theme and workstream colours are independent. Combine them in one URL, such as `?theme=sand&palette=muted`.

The chart bands connect successive proposed gate dates. The first band starts at the chart's
September boundary for display. These bands do not specify task start dates or effort.
The tracker dates remain in the snapshot for reference. The chart uses the outline dates.

Task counts describe the linked records. They do not establish gate completion.
The source notes identify differences between the outline and the tracker.

## Refresh

The Project dashboard workflow refreshes the tracker fields every six hours and on manual dispatch
from `main`. `fetch-snapshot.mjs` queries the tracker data source through the Notion API and
replaces the epics and tasks in `snapshot.json` (titles, statuses, owners, dates, URLs) on the
runner. The workflow then builds and deploys from that refreshed copy. Nothing is committed: the
committed snapshot is the hand-edited baseline, and the published page shows the tracker as of the
last run. A failed API call or a snapshot that fails validation stops the run, and the published
page keeps its last deployment. The iterations, the workstream names and colours, and the proposed
iteration mapping stay as the snapshot holds them. A new task takes the first iteration of its epic.

The workflow needs the repository secret `NOTION_TOKEN`: a Notion token that can read the tracker
database. To refresh the committed baseline locally:

```sh
NOTION_TOKEN=... node dashboard/fetch-snapshot.mjs && node dashboard/build.mjs
```

Iteration statuses, notes, and dates come from the presentation outline. Edit them by hand in
`snapshot.json` when the outline changes, and set `outlineEditedAt` to the outline's edit time.

The Notion attachment is a separate snapshot. To update it, upload the generated HTML through
the Notion connector and replace the page's embed with the returned `file-upload://` source.

## GitHub Pages

The Project dashboard workflow builds the site and publishes `index.html` and the query sequence
diagram, `docs/diagrams/census-query-sequence.html`, as a second page. The diagram is an Archify
viewer with a light and dark theme. It checks that Pages has private access before each deployment. GitHub sign-in and repository
read access are required. The internal repository grants read access to enterprise members.

The workflow deploys dashboard changes from `feat/project-dashboard` and `main`.
The feature branch supports the initial deployment before merge. Delete that branch after merge
to use `main` for subsequent updates. Pull requests run the build check without deployment.

The `github-pages` environment permits only those two branches. The Pages source is GitHub
Actions, and the Pages API must report `public: false`.
Repository visibility alone does not establish the site's access controls.
See [GitHub Pages access controls](https://docs.github.com/en/enterprise-cloud@latest/pages/getting-started-with-github-pages/changing-the-visibility-of-your-github-pages-site).
