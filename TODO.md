# Performance Follow-Up

- The `dbh` slider lag improved after debouncing and removing unrelated landscape refreshes.
- The remaining noticeable lag is the first visit to `Treatment Report`.
- Measured treatment execution was not the main bottleneck. `do_treatment()` was roughly `80-90 ms`.
- `Cut Report` was acceptable at roughly `30-45 ms`.
- The dominant remaining cost was `StandViewCoordinator.update_treatment_report`, measured around `800-890 ms`.
- Within that path, clump size computation was a meaningful but not dominant contributor:
  - `get_clump_sizes()` was roughly `120-130 ms`
  - `get_treat_clump_sizes()` was roughly `45 ms` after treatment
- Most of the remaining cost is likely in Python-side treatment report work:
  - repeated TAO-to-NumPy conversion and basal area prep
  - four density plot rebuilds
  - `gaussian_kde`
  - drawing four Matplotlib canvases

## Next Optimization Targets

- Avoid full treatment-report plot recomputation on target edits when only labels or text changed.
- Add finer-grained profiling inside `update_treatment_report()` if more measurements are needed.
- Reuse current-stand distributions within the report lifecycle instead of recomputing them unnecessarily.
- Consider replacing KDE-based clump-size plots with a cheaper histogram or smoothed binned curve if UI responsiveness remains poor.
