# md_box_trend.awk -- is the bilayer still changing size?
#
# Reads a gmx energy xvg with Box-X and Box-Y columns and compares the mean
# xy area over the first quarter of the trajectory with the last quarter.
#
# Deliberately NOT reported: area per lipid. Dividing the box xy area by the
# lipid count attributes the protein's own cross-section to the lipids, which
# inflates the figure past any experimental POPC value. It would look like a
# meaningful number and would not be one. What the box honestly answers is
# whether the membrane has stopped changing size: first and last quarters that
# agree mean equilibrated; a still-shrinking box means it is not.

/^[@#]/ { next }
{ n++; a[n] = $2 * $3 }
END {
    if (n < 8) {
        printf "    bilayer  too few frames to judge drift\n"
        exit
    }
    q = int(n / 4)
    for (i = 1; i <= q; i++) first += a[i]
    for (i = n - q + 1; i <= n; i++) last += a[i]
    first /= q
    last /= q
    printf "    bilayer  box xy area %.2f -> %.2f nm^2 (first vs last quarter), %+.2f%%\n", first, last, 100 * (last - first) / first
    printf "             area per lipid not reported: the protein occupies part of this area\n"
}
