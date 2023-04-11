# OnishRDB

Osborne Nishimura R Database interface

## Requirements

Libraries R.cache, RMariaDb, dplyr.

``` r
if (! require(R.cache)) {
  install.packages('R.cache')
  library(R.cache)
}
if (! require(RMariaDB)) {
  install.packages('RMariaDB')
  library(RMariaDB)
}
if (! require(dplyr)) {
  install.packages('dplyr')
  library(dplyr)
}
```
