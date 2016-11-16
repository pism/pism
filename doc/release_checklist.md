### PISM release checklist

1.  [ ] Create a new `stableXX` branch, if needed.
2.  [ ] Merge, if needed.
3.  [ ] Set `Pism_BRANCH` in `CMakeLists.txt`
4.  [ ] Replace `0.x` with `0.x+1` in docs.
7.  [ ] Update links to `INSTALL.md` in docs.
5.  [ ] Update `CHANGES.md`.
6.  [ ] Tag.
```
git tag -a v0.X -m "The v0.X release. See CHANGES.md for the list of changes since v0.X-1."
```
7.  [ ] Push.
```
git push -u origin HEAD
```
8.  [ ] Push tags.
```
git push --tags
```
9.  [ ] Re-build docs.
```
make manual installation forcing browser.tgz pismpython_docs
```
10. [ ] Upload these docs.
11. [ ] Switch the default branch on [github.com/pism/pism](http:github.com/pism/pism)
12. [ ] Write a news item for `pism-docs.org`.
13. [ ] Replace `0.x` with `0.x+1` on `pism-docs.org`.
14. [ ] Send an e-mail to `pismusers@pismdocs.org`.
15. [ ] Tell more people, if desired.
