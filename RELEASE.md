# Release Process

This project uses Semantic Versioning: `MAJOR.MINOR.PATCH`.

- Increment `PATCH` for bug fixes and documentation-only corrections.
- Increment `MINOR` for backward-compatible features.
- Increment `MAJOR` for breaking API or behavior changes.

## Checklist

1. Choose the next version number.
2. Update `__version__` in `cancermuts/__init__.py`.
3. Move entries from `CHANGELOG.md` under `[Unreleased]` into a new dated release section.
4. Commit the release changes:

   ```bash
   git add cancermuts/__init__.py setup.py CHANGELOG.md RELEASE.md
   git commit -m "Release X.Y.Z"
   ```

5. Create an annotated tag:

   ```bash
   git tag -a vX.Y.Z -m "Release X.Y.Z"
   ```

6. Push the branch and tag:

   ```bash
   git push
   git push origin vX.Y.Z
   ```

GitHub releases can then be created from the version tag, using the matching `CHANGELOG.md` section as the release notes.
