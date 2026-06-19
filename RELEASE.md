# Release Process

This project uses Semantic Versioning: `MAJOR.MINOR.PATCH`.

- Increment `PATCH` for bug fixes and documentation-only corrections.
- Increment `MINOR` for backward-compatible features.
- Increment `MAJOR` for breaking API or behavior changes.

The `master` branch should always represent the latest released version. Unreleased work is collected on the long-lived `devel` branch until it is ready to become the next release.

## Branch Workflow

1. Branch feature and fix work from `devel`.
2. Open a pull request from each completed feature or fix branch into `devel`.
3. Review and merge the pull request into `devel`.
4. Record release notes in `CHANGELOG.md` under `[Unreleased]` while work is still in progress.
5. Keep `__version__` in `cancermuts/__init__.py` at the latest released version while `devel` is accumulating unreleased changes.
6. When `devel` is ready to release, choose the next version number, for example `2.1.0` or `2.0.1`.
7. Update `__version__` in `cancermuts/__init__.py`.
8. Move the relevant `CHANGELOG.md` entries from `[Unreleased]` into a new dated release section:

   ```markdown
   ## [X.Y.Z] - YYYY-MM-DD
   ```

9. Commit the release preparation changes on `devel`:

   ```bash
   git add cancermuts/__init__.py setup.py CHANGELOG.md RELEASE.md
   git commit -m "Release X.Y.Z"
   ```

10. Merge `devel` into `master`:

   ```bash
   git switch master
   git merge devel
   ```

11. Create an annotated tag on the released `master` commit:

   ```bash
   git tag -a vX.Y.Z -m "Release X.Y.Z"
   ```

12. Push `master`, `devel`, and the tag:

   ```bash
   git push origin master devel
   git push origin vX.Y.Z
   ```

GitHub releases can then be created from the version tag, using the matching `CHANGELOG.md` section as the release notes.

## Hotfixes

For a patch release that cannot wait for the rest of `devel`, branch from `master`, apply only the fix, update the version and changelog, then merge and tag it in the same way:

```bash
git switch -c hotfix_X.Y.Z master
```

After the hotfix is released, merge `master` back into `devel` so the rolling branch includes the patch release.
