# Census 2 project dashboard

The dashboard shows the iteration gates from the presentation outline. Each iteration links to
the relevant Notion epics and tasks. The dashboard uses a dated snapshot of the tracker.

[GitHub Pages hosts the team dashboard](https://scaling-dollop-pzyy6pk.pages.github.io/)
with private access. The
[Notion dashboard page](https://app.notion.com/p/3ea16a75378981019542e9e2a3db2150)
also holds the original HTML snapshot.

Open `index.html` in a browser. The file also works as an HTML attachment in Notion.
It contains its data, styles, and script. It makes no API requests and needs no credentials.

## Build

Run this command from the repository root. It needs Node.js 18 or later.

```sh
node dashboard/build.mjs
```

The build validates the snapshot and writes `dashboard/index.html`. Commit the generated file
with the source files so the preview opens without a build step.

## Sources and dates

- The [presentation outline](https://app.notion.com/p/3e616a753789813491d6ffa66a4c3998)
  supplies the iteration names, proposed gate dates, and outline statuses.
- The [project](https://app.notion.com/p/38f16a753789806b8fb0e21950fd3659)
  supplies the epic membership.
- The [tracker](https://app.notion.com/p/95716a7537898228b0150198eacd21dc)
  supplies task titles, statuses, owners, and parent relations.
- The [plan](https://app.notion.com/p/3b216a75378981c9b7f3d954789aefe4)
  and [iteration proposal](https://app.notion.com/p/68916a75378982a3908e018e928d7eb7)
  supply the phase deliverables used for the proposed iteration mapping.

The tracker has no iteration property. `snapshot.json` holds the proposed mapping.
An epic can support several iterations. Each child task has one proposed iteration.
The dashboard labels this mapping for review. It does not change the tracker.

Each epic has an original colour in `snapshot.json`. Phase numbers identify the chart bands.
The bands do not assign separate epic dates.

Midnight is the default theme. A saved choice or URL theme takes precedence.
The theme and epic colour pickers start closed. Click the colour-settings icon to open them.
The dashboard theme selector offers Slate, Sand, and Midnight. System uses the original green theme
and follows the device's light or dark setting. The browser saves the choice when storage is
available. A URL with `?theme=slate`, `?theme=sand`, or `?theme=midnight` opens that theme.

The epic colour selector compares the current palette with muted and vivid variants from `palettes.json`.
The selection applies to the chart, legend, and links. The browser saves the choice when storage
is available. A URL with `?palette=muted` or `?palette=vivid` opens the selected variant.
Theme and epic colours are independent. Combine them in one URL, such as `?theme=sand&palette=muted`.

The chart bands connect successive proposed gate dates. The first band starts at the chart's
September boundary for display. These bands do not specify task start dates or effort.
The tracker dates remain in the snapshot for reference. The chart uses the outline dates.

Task counts describe the linked records. They do not establish gate completion.
The source notes identify differences between the outline and the tracker.

## Refresh

The published HTML does not refresh automatically.

1. Fetch the outline and the project page through the Notion connector.
2. Fetch the tracker schema for data source `b0816a75-3789-82f1-91e3-87989b8c7954`.
3. Query the project's epics and their child tasks by the parent relation.
4. Update the source fields in `snapshot.json` and its `capturedAt` timestamp.
5. Review the proposed iteration mapping and the source notes.
6. Run the build command.
7. Open the preview and check the iteration selection, task filters, and links.
8. Commit the updated source and generated HTML.
9. Push the deployment branch and check the Project dashboard workflow.

The Notion attachment is a separate snapshot. To update it, upload the generated HTML through
the Notion connector and replace the page's embed with the returned `file-upload://` source.

## GitHub Pages

The Project dashboard workflow builds the site and publishes only `index.html`.
It checks that Pages has private access before each deployment. GitHub sign-in and repository
read access are required. The internal repository grants read access to enterprise members.

The workflow deploys dashboard changes from `feat/project-dashboard` and `main`.
The feature branch supports the initial deployment before merge. Delete that branch after merge
to use `main` for subsequent updates. Pull requests run the build check without deployment.

The `github-pages` environment permits only those two branches. The Pages source is GitHub
Actions, and the Pages API must report `public: false`.
Repository visibility alone does not establish the site's access controls.
See [GitHub Pages access controls](https://docs.github.com/en/enterprise-cloud@latest/pages/getting-started-with-github-pages/changing-the-visibility-of-your-github-pages-site).
