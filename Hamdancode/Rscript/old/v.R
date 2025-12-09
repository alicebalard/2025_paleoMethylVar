# 1) See where R is looking for packages
.libPaths()

# 2) Remove any leftover lock dirs (sometimes cause corruption on Windows)
unlink(Sys.glob(file.path(.libPaths()[1], "00LOCK*")), recursive = TRUE, force = TRUE)

# 3) Remove gtable from *all* library paths, just in case
sapply(.libPaths(), function(p) unlink(file.path(p, "gtable"), recursive = TRUE, force = TRUE))
suppressWarnings(remove.packages("gtable"))

# 4) Reinstall a fresh binary build (fastest on Windows)
install.packages("gtable", type = "binary")

# 5) Verify it loads
library(gtable)
packageVersion("gtable")
