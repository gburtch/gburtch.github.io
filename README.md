# gburtch.github.io

Personal site for Gord Burtch, built on Jekyll with the
[Academic Pages](https://github.com/academicpages/academicpages.github.io)
fork of Minimal Mistakes (MIT, see LICENSE). Hosted on GitHub Pages from `master`.

## Layout

| Path | What it is |
| --- | --- |
| `_pages/about.md` | Home page / bio (`/`) |
| `_pages/cv.md` | CV page; links to the PDF below |
| `files/Gordon_Burtch_CV/` | CV source (`gb_cv.tex`) and compiled `gb_cv.pdf` |
| `_posts/` | Blog posts, `YYYY-MM-DD-slug.md` |
| `_drafts/` | R Markdown sources for posts |
| `images/<post>_files/` | Figures referenced by posts |
| `_data/navigation.yml` | Top navigation |
| `_config.yml` | Site settings |

## Adding a post

1. Write the post in R Markdown under `_drafts/` and knit to GitHub-flavoured
   markdown (`output: github_document`), or write markdown directly.
2. Save it as `_posts/YYYY-MM-DD-slug.md` with front matter:

   ```yaml
   ---
   title: "Post title"
   date: "YYYY-MM-DD"
   permalink: /posts/YYYY/MM/slug/
   tags:
     - econometrics
   ---
   ```

3. Put figures in `images/<slug>_files/` and reference them as
   `/images/<slug>_files/...`. MathJax is enabled, so `$...$` works.
4. Commit and push to `master`; GitHub Pages rebuilds in a minute or two.

## Running locally

```bash
bundle install
bundle exec jekyll serve
```

Then open http://localhost:4000.
