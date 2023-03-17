Lab 10: Advanced SSF Models - integrated SSF’s
================
Mark Hebblewhite and Eric Palm
March 17, 2023

## Preliminaries: setting packages

``` r
#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


packages <- c("sf","terra","lubridate", "tidyverse","ggplot2","mapview","maptools","leaflet","xtable","broom","stars","magrittr","cowplot", "tmap","suncalc", "survival", "amt", "glmmTMB", "TMB")

#run function to install packages
ipak(packages)
```

A note about tidyverse. Tidyverse is REALLY useful, but, sometimes has
‘hidden’ conflicts with other packages. Check which specific conflicts
there are with:

``` r
tidyverse_conflicts()
```

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::as.difftime() masks base::as.difftime()
    ## ✖ lubridate::date()        masks base::date()
    ## ✖ magrittr::extract()      masks tidyr::extract(), terra::extract()
    ## ✖ amt::filter()            masks dplyr::filter(), stats::filter()
    ## ✖ magrittr::inset()        masks terra::inset()
    ## ✖ lubridate::intersect()   masks terra::intersect(), base::intersect()
    ## ✖ dplyr::lag()             masks stats::lag()
    ## ✖ magrittr::set_names()    masks purrr::set_names()
    ## ✖ lubridate::setdiff()     masks base::setdiff()
    ## ✖ cowplot::stamp()         masks lubridate::stamp()
    ## ✖ lubridate::union()       masks terra::union(), base::union()

Note especially this line:

    amt::select()            masks raster::select(), dplyr::select()

Thus, for any tidyverse nesting and selecting, we need to specify
dplyr::select().

# Lab 10 - Advanced Step Selection Function Models

In today lab we will continue to build skills modeling resource
selection using movement-based statistical models. In particular, we
will continue with the amt package for a case study of Fishers, and fit
integrated step selection functions (iSSF) to individual and multiple
individuals from a well known Fisher dataset. Today’s lab objectives
are:

1)  Fit an iSSF to a single Fisher dataset, estimating a clogit model,
    as well as then fitting the resultant movement and habitat kernels
    to field data.
2)  Fit an iSSF to multiple individual Fishers.
3)  Explore mixed-effects cLogit models using the coxme and mcclogit
    packages.
4)  Illustrate the recommended way to fit clogit models using glmmTMB
    based on Muff et al. (2020)

# Animal movement tools (amt)

Next, we will start using the amt:: R package for managing tracking data
and conducting habitat selection analyses.

We will estimate a Step Selection Function, and, an integrated SSF to
map the predictions of spatial distribution from the SSF based on the
work of Tal Avgar et al’s integrated Step Selection Function modeling
approach. This will be based on the examples provided in this paper by:

Singer, J., Fieberg, J. & Avgar, T. (2019) Animal movement tools (amt):
R package for managing tracking data and conducting habitat selection
analyses. Ecol Evol, 9, 880-890.https://doi.org/10.1002/ece3.4823

based on this method:

Avgar, T., Potts, J.R., Lewis, M.A., Boyce, M.S. & Börger, L. (2016)
Integrated step selection analysis: bridging the gap between resource
selection and animal movement. Methods in Ecology and Evolution, 7,
619-630. <https://doi.org/10.1111/2041-210X.12528>

These data are from Fishers in upstate New York collected by a
colleague, Scott La Point, during his PhD with Max Planck.  
![Figure 10.1 Fisher](Figures/fisher.jpg)

To get started, we will review the handy Figure from Avgar et al. about
the steps involved in fitting an integrated Step Selection Function from
animal movement data.

![Figure 10.2. iSSF workflow](Figures/AvgarFig1.jpg)

Lets get started….

``` r
fisher <- read.csv("Data/Martes pennanti LaPoint New York.csv")
head(fisher)
```

    ##    event.id visible               timestamp location.long location.lat
    ## 1 170635808    true 2011-02-11 17:30:00.000            NA           NA
    ## 2 170635809    true 2011-02-11 17:40:00.000            NA           NA
    ## 3 170635810    true 2011-02-11 17:50:00.000            NA           NA
    ## 4 170635811    true 2011-02-11 18:00:00.000            NA           NA
    ## 5 170635812    true 2011-02-11 18:06:13.999     -73.86015     42.79528
    ## 6 170635813    true 2011-02-11 18:10:08.999     -73.86001     42.79534
    ##   behavioural.classification eobs.battery.voltage eobs.fix.battery.voltage
    ## 1                                              NA                       NA
    ## 2                                              NA                       NA
    ## 3                                              NA                       NA
    ## 4                                              NA                       NA
    ## 5                                              NA                       NA
    ## 6                                              NA                       NA
    ##   eobs.horizontal.accuracy.estimate eobs.key.bin.checksum
    ## 1                                NA                    NA
    ## 2                                NA                    NA
    ## 3                                NA                    NA
    ## 4                                NA                    NA
    ## 5                                NA                    NA
    ## 6                                NA                    NA
    ##   eobs.speed.accuracy.estimate eobs.start.timestamp eobs.status
    ## 1                           NA                   NA          NA
    ## 2                           NA                   NA          NA
    ## 3                           NA                   NA          NA
    ## 4                           NA                   NA          NA
    ## 5                           NA                   NA          NA
    ## 6                           NA                   NA          NA
    ##   eobs.temperature eobs.type.of.fix eobs.used.time.to.get.fix ground.speed
    ## 1               NA               NA                       110           NA
    ## 2               NA               NA                       110           NA
    ## 3               NA               NA                       110           NA
    ## 4               NA               NA                       110           NA
    ## 5               NA               NA                        80           NA
    ## 6               NA               NA                         7           NA
    ##   heading height.above.ellipsoid manually.marked.outlier sensor.type
    ## 1      NA                     NA                      NA         gps
    ## 2      NA                     NA                      NA         gps
    ## 3      NA                     NA                      NA         gps
    ## 4      NA                     NA                      NA         gps
    ## 5      NA                     NA                      NA         gps
    ## 6      NA                     NA                      NA         gps
    ##   individual.taxon.canonical.name tag.local.identifier
    ## 1                 Martes pennanti                 1072
    ## 2                 Martes pennanti                 1072
    ## 3                 Martes pennanti                 1072
    ## 4                 Martes pennanti                 1072
    ## 5                 Martes pennanti                 1072
    ## 6                 Martes pennanti                 1072
    ##   individual.local.identifier                       study.name utm.easting
    ## 1                          F1 Martes pennanti LaPoint New York          NA
    ## 2                          F1 Martes pennanti LaPoint New York          NA
    ## 3                          F1 Martes pennanti LaPoint New York          NA
    ## 4                          F1 Martes pennanti LaPoint New York          NA
    ## 5                          F1 Martes pennanti LaPoint New York    593215.5
    ## 6                          F1 Martes pennanti LaPoint New York    593226.9
    ##   utm.northing utm.zone        study.timezone   study.local.timestamp
    ## 1           NA          Eastern Standard Time 2011-02-11 12:30:00.000
    ## 2           NA          Eastern Standard Time 2011-02-11 12:40:00.000
    ## 3           NA          Eastern Standard Time 2011-02-11 12:50:00.000
    ## 4           NA          Eastern Standard Time 2011-02-11 13:00:00.000
    ## 5      4738712      18N Eastern Standard Time 2011-02-11 13:06:13.999
    ## 6      4738718      18N Eastern Standard Time 2011-02-11 13:10:08.999

``` r
str(fisher)
```

    ## 'data.frame':    47347 obs. of  30 variables:
    ##  $ event.id                         : int  170635808 170635809 170635810 170635811 170635812 170635813 170635814 170635815 170635816 170635817 ...
    ##  $ visible                          : chr  "true" "true" "true" "true" ...
    ##  $ timestamp                        : chr  "2011-02-11 17:30:00.000" "2011-02-11 17:40:00.000" "2011-02-11 17:50:00.000" "2011-02-11 18:00:00.000" ...
    ##  $ location.long                    : num  NA NA NA NA -73.9 ...
    ##  $ location.lat                     : num  NA NA NA NA 42.8 ...
    ##  $ behavioural.classification       : chr  "" "" "" "" ...
    ##  $ eobs.battery.voltage             : logi  NA NA NA NA NA NA ...
    ##  $ eobs.fix.battery.voltage         : logi  NA NA NA NA NA NA ...
    ##  $ eobs.horizontal.accuracy.estimate: logi  NA NA NA NA NA NA ...
    ##  $ eobs.key.bin.checksum            : logi  NA NA NA NA NA NA ...
    ##  $ eobs.speed.accuracy.estimate     : logi  NA NA NA NA NA NA ...
    ##  $ eobs.start.timestamp             : logi  NA NA NA NA NA NA ...
    ##  $ eobs.status                      : logi  NA NA NA NA NA NA ...
    ##  $ eobs.temperature                 : logi  NA NA NA NA NA NA ...
    ##  $ eobs.type.of.fix                 : logi  NA NA NA NA NA NA ...
    ##  $ eobs.used.time.to.get.fix        : int  110 110 110 110 80 7 58 110 90 110 ...
    ##  $ ground.speed                     : logi  NA NA NA NA NA NA ...
    ##  $ heading                          : logi  NA NA NA NA NA NA ...
    ##  $ height.above.ellipsoid           : logi  NA NA NA NA NA NA ...
    ##  $ manually.marked.outlier          : logi  NA NA NA NA NA NA ...
    ##  $ sensor.type                      : chr  "gps" "gps" "gps" "gps" ...
    ##  $ individual.taxon.canonical.name  : chr  "Martes pennanti" "Martes pennanti" "Martes pennanti" "Martes pennanti" ...
    ##  $ tag.local.identifier             : int  1072 1072 1072 1072 1072 1072 1072 1072 1072 1072 ...
    ##  $ individual.local.identifier      : chr  "F1" "F1" "F1" "F1" ...
    ##  $ study.name                       : chr  "Martes pennanti LaPoint New York" "Martes pennanti LaPoint New York" "Martes pennanti LaPoint New York" "Martes pennanti LaPoint New York" ...
    ##  $ utm.easting                      : num  NA NA NA NA 593215 ...
    ##  $ utm.northing                     : num  NA NA NA NA 4738712 ...
    ##  $ utm.zone                         : chr  "" "" "" "" ...
    ##  $ study.timezone                   : chr  "Eastern Standard Time" "Eastern Standard Time" "Eastern Standard Time" "Eastern Standard Time" ...
    ##  $ study.local.timestamp            : chr  "2011-02-11 12:30:00.000" "2011-02-11 12:40:00.000" "2011-02-11 12:50:00.000" "2011-02-11 13:00:00.000" ...

``` r
ggplot(fisher, aes(location.long, location.lat, colour = individual.local.identifier)) + geom_point()
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Next, let’s turn fisher into a sf data frame for use in raster
operations later. First we need to remove NA’s, which we can do here
using complete.cases and which we do later using the filter(!is.na())
command.

``` r
fisher1<-fisher[complete.cases(fisher[4:5]),]
fisher_sf <- st_as_sf(fisher1,
                      coords = c("location.long","location.lat"),
                      crs = "EPSG:4326")
```

``` r
mapview(fisher_sf, zcol="individual.local.identifier", legend = TRUE, cex=5, lwd=2, map.type = c("OpenStreetMap.DE", "Esri.WorldShadedRelief"))
```

![Fisher Map](Figures/Rplot_line95.png) What is always amazing to me is
how urban these Fisher’s are. Whereas in the west, Fishers tend to be
more of a wilderness species. These are living in downtown urban center
parks making a tidy living eating gray squirrels, I bet. Next we will
bringing Fisher data from a single individual, ID 1016, Fisher M1 (known
as “RickyT”), into a MOVE object.

``` r
dat <- read_csv("Data/Martes pennanti LaPoint New York.csv") %>%
   filter(!is.na(`location-lat`)) %>%
   dplyr::select(x = `location-long`, y = `location-lat`,
           t = `timestamp`, id = `tag-local-identifier`) %>%
    filter(id %in% c(1465, 1466, 1072, 1078, 1016, 1469)) # for example 2
dat_1 <- dat %>% filter(id == 1016)
head(dat_1)
```

    ## # A tibble: 6 × 4
    ##       x     y t                      id
    ##   <dbl> <dbl> <dttm>              <dbl>
    ## 1 -73.9  42.8 2010-02-09 17:01:23  1016
    ## 2 -73.9  42.8 2010-02-09 17:28:56  1016
    ## 3 -73.9  42.8 2010-02-09 17:31:28  1016
    ## 4 -73.9  42.8 2010-02-09 17:33:22  1016
    ## 5 -73.9  42.8 2010-02-09 17:35:34  1016
    ## 6 -73.9  42.8 2010-02-09 17:36:10  1016

``` r
ggplot(dat_1, aes(x, y, colour = id)) + geom_point()
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

*From Singer et al. (2019)…* The function amt::make_track creates a
track (the basic building block of the amt package), given the names of
the columns containing x and y coordinates, time (t), and we can set a
coordinate reference system (CRS). The original data were provided in
geographical coordinates (EPSG code: 4326). Here, we shall transform
this original CRS (using function amt::transform \_ coords).

We will make a track just with the data from RickyT.

``` r
dat_R <- amt::make_track(dat_1, x, y, t, crs = 4326) %>%
     amt::transform_coords(3857)


amt::summarize_sampling_rate(dat_R)
```

    ## # A tibble: 1 × 9
    ##     min    q1 median  mean    q3   max    sd     n unit 
    ##   <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr>
    ## 1 0.100  1.93   2.03  8.04  2.57 1209.  44.0  8957 min

*From Singer et al. (2019)..*

We see that we have 8,957 total locations, the shortest interval between
locations is 0.1 min and the largest time interval between locations is
1,208 min, with median interval length equal to roughly 2 min. Despite
the 2 min temporal resolution, we choose to resample the track to 10 min
with a tolerance of 1 min (amt::track_resample), to conduct the analyses
on the same temporal scale as the next example (where most individuals
had a median sampling rate of 10 min).

The function minutes from the package tidyverse and lubridate (Grolemund
& Wickham, 2011), is used here to create an object of class Period that
is then passed to amt::track_resample. Periods can be specified using
all common time units; thus, it is straightforward to specify a sampling
rate and an acceptable tolerance. We will also choose to keep only those
bursts (subsets of the track with constant sampling rate, within the
specified tolerance) with at least three relocations, the minimum
required to calculate a turn angle amt::filter_min_n\_burst). Note a
burst is one set of consistent locations without a missed fix.

The following code implements those choices and translates from a point
representation to a step (step length, turn angle) representation of the
data. In the final line of the code snippet, we use the function
amt::time_of_day (a wrapper around maptools:: sunriset and
maptools::crepuscule; Bivand & Lewin‐Koh, 2017) to calculate if a
location was taken during the day or night. If the argument
include.crepuscule is set to TRUE, the function not only considers day
and night, but also dawn and dusk.

``` r
stps <- amt::track_resample(dat_R, rate = minutes(10), tolerance = minutes(1)) %>%
  amt::filter_min_n_burst(min_n = 3) %>% amt::steps_by_burst() %>%
  amt::time_of_day(include.crepuscule = FALSE)

str(stps, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)
```

    ## brstd_s_ [1,494 × 12] (S3: bursted_steps_xyt/steps_xyt/steps_xy/tbl_df/tbl/data.frame)
    ##  $ burst_     : num [1:1494] 16 16 16 16 18 18 21 21 23 23 ...
    ##  $ x1_        : Named num [1:1494] -8226199 -8226180 -8226179 -8226170 -8226177 ...
    ##  $ x2_        : Named num [1:1494] -8226180 -8226179 -8226170 -8226182 -8226182 ...
    ##  $ y1_        : Named num [1:1494] 5287115 5286844 5286852 5286826 5286843 ...
    ##  $ y2_        : Named num [1:1494] 5286844 5286852 5286826 5286835 5286835 ...
    ##  $ sl_        : Named num [1:1494] 271.01 8.17 28.29 15.23 9.33 ...
    ##  $ direction_p: Named num [1:1494] -1.5 1.46 -1.23 2.49 -2.11 ...
    ##  $ ta_        : Named num [1:1494] NA 2.96 -2.7 -2.56 NA ...
    ##  $ t1_        : POSIXct[1:1494], format: "2010-02-10 03:00:35" "2010-02-10 03:10:39" ...
    ##  $ t2_        : POSIXct[1:1494], format: "2010-02-10 03:10:39" "2010-02-10 03:20:10" ...
    ##  $ dt_        : 'difftime' num [1:1494] 10.0666666666667 9.50001666545868 10.2333333333333 10.3666833321253 ...
    ##  $ tod_end_   : Factor w/ 2 levels "day","night": 2 2 2 2 2 2 2 2 2 2 ...

## Obtaining the NLCD data from `FedData`

Next, we will obtain the Landcover data from the National Landcover
Database using the FedData pacakge. I encourage you to explore the
FedData package here - they include several very handy spatial datasets
including NLCD, Daymet, and national hydrography datasets. This is a
VERY handy R package.

See here <https://cran.r-project.org/web/packages/FedData/FedData.pdf>
and <https://github.com/ropensci/FedData> and Kyle’s website here:
<https://www.bocinsky.io>

The first set of commands are to create a blank raster that we can use
as an extent to ‘clip’ the incoming NLCD data by. We also need to make
the projection a PROJECTED, cartesian UTM NAD83 projection as that is
what the NLCD are stored in.

First, we have to install the R package `FedData` from github because of
changes to NLCD file services which are now being delivered using a
different web coverage service.

``` r
require("devtools")
devtools::install_github("ropensci/FedData")
library(FedData)
```

The key here is to remind yourself that the extent will have to be
bigger than the actual extent of the Fisher data, because of the
generation of availability samples during the SSF models. Keep this in
mind in your own research as well. For example, I had to go all the way
back to the beginning and clip out a larger extent later.

``` r
## Create a Mask Raster based on the extent of Fisher. Note that I made it arbitrarily larger. 
crs(fisher_sf, proj=TRUE)
```

    ## [1] "+proj=longlat +datum=WGS84 +no_defs"

``` r
fisher_sf2 <-st_transform(fisher_sf, crs = "EPSG:3857")
ext(fisher_sf2)
```

    ## SpatExtent : -8230998.57111662, -8169678.42998823, 5267780.0505195, 5289319.21290995 (xmin, xmax, ymin, ymax)

``` r
st_bbox(fisher_sf2)
```

    ##     xmin     ymin     xmax     ymax 
    ## -8230999  5267780 -8169678  5289319

``` r
mask.raster<-rast()
ext(mask.raster) <- c(xmin=-8240000, xmax=-8160000, ymin=5260000 , ymax=5295000)

res(mask.raster) = 30
crs(mask.raster)<- "+init=epsg:3857"  
#set all values of mask.raster to zero
mask.raster[]<-0
#projection(mask.raster)
plot(mask.raster)
plot(fisher_sf2, add = TRUE)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- --> This
confirms that our mask contains all the Fisher data, and, a buffer that
is probably (?) large enough to contain all AVAILABLE samples created by
the ssf available point generation method. This is a common rookie
mistake. Always remember to buffer your used points.

Next, we use the get_nlcd command to obtain 1 tile of the NLCD dataset
based on our input mask.raster. The full command can be specified with
options such as the location for the raw.dir, the extraction.dir,
raster.options, and whether you want to ‘overwrite’(force redo). For
example, here:

    get_nlcd(mask.raster, label="landuse", year = 2011, dataset = "landcover", raw.dir = "/Users/mark.hebblewhite/Box Sync/Teaching/UofMcourses/WILD562/Spring2021/Labs/Lab10/Data/raw/", extraction.dir = "/Users/mark.hebblewhite/Box Sync/Teaching/UofMcourses/WILD562/Spring2021/Labs/Lab10/2021/Fisher/extract/", raster.options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), force.redo = F)

Note, here, I do create a directory on my local machine where to put the
raw downloaded NLCD tiles and extractions.

``` r
get_nlcd(mask.raster, label = "landuse", year = 2016, dataset = "landcover", extraction.dir = "Data/extract/", force.redo = TRUE)
```

    ## class      : RasterLayer 
    ## dimensions : 1290, 2096, 2703840  (nrow, ncol, ncell)
    ## resolution : 30, 30  (x, y)
    ## extent     : 1769265, 1832145, 2390655, 2429355  (xmin, xmax, ymin, ymax)
    ## crs        : +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
    ## source     : landuse_NLCD_Land_Cover_2016.tif 
    ## names      : Class 
    ## values     : 11, 95  (min, max)
    ## attributes :
    ##         ID Class   Color Description
    ##  from:   0  <NA> #000000        <NA>
    ##   to : 255  <NA> #000000        <NA>

``` r
land_use <- rast("Data/extract/landuse_NLCD_Land_Cover_2016.tif")
str(values(land_use))
```

    ##  int [1:2703840, 1] 90 90 90 90 90 90 90 90 41 41 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : NULL
    ##   ..$ : chr "Class"

``` r
head(land_use)
```

    ##            Class
    ## 1 Woody Wetlands
    ## 2 Woody Wetlands
    ## 3 Woody Wetlands
    ## 4 Woody Wetlands
    ## 5 Woody Wetlands
    ## 6 Woody Wetlands

``` r
hist(values(land_use))
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Like any landcover model, there are a number of categories, \~ 250 for
the Conterminous United States (CONUS). In our upstate New York study
area, we have 17.

``` r
unique(values(land_use))
```

    ##       Class
    ##  [1,]    90
    ##  [2,]    41
    ##  [3,]    43
    ##  [4,]    81
    ##  [5,]    95
    ##  [6,]    23
    ##  [7,]    22
    ##  [8,]    21
    ##  [9,]    82
    ## [10,]    42
    ## [11,]    71
    ## [12,]    24
    ## [13,]    11
    ## [14,]    52
    ## [15,]    31

``` r
rat <- as.data.frame(cats(land_use)) %>% 
  filter(!is.na(Class))
rat
```

    ##    ID                        Class   Color
    ## 1  11                   Open Water #5475A8
    ## 2  12           Perennial Ice/Snow #FFFFFF
    ## 3  21        Developed, Open Space #E8D1D1
    ## 4  22     Developed, Low Intensity #E29E8C
    ## 5  23  Developed, Medium Intensity #ff0000
    ## 6  24     Developed High Intensity #B50000
    ## 7  31 Barren Land (Rock/Sand/Clay) #D2CDC0
    ## 8  41             Deciduous Forest #85C77E
    ## 9  42             Evergreen Forest #38814E
    ## 10 43                 Mixed Forest #D4E7B0
    ## 11 51                  Dwarf Scrub #AF963C
    ## 12 52                  Shrub/Scrub #DCCA8F
    ## 13 71         Grassland/Herbaceous #FDE9AA
    ## 14 72             Sedge/Herbaceous #D1D182
    ## 15 73                      Lichens #A3CC51
    ## 16 74                         Moss #82BA9E
    ## 17 81                  Pasture/Hay #FBF65D
    ## 18 82             Cultivated Crops #CA9146
    ## 19 90               Woody Wetlands #C8E6F8
    ## 20 95 Emergent Herbaceous Wetlands #64B3D5
    ##                                                                                                                                                                                                                                                                                                                                                                Description
    ## 1                                                                                                                                                                                                                                                                                           Areas of open water, generally with less than 25% cover of vegetation or soil.
    ## 2                                                                                                                                                                                                                                                                  Areas characterized by a perennial cover of ice and/or snow, generally greater than 25% of total cover.
    ## 3  Areas with a mixture of some constructed materials, but mostly vegetation in the form of lawn grasses. Impervious surfaces account for less than 20% of total cover. These areas most commonly include large-lot single-family housing units, parks, golf courses, and vegetation planted in developed settings for recreation, erosion control, or aesthetic purposes.
    ## 4                                                                                                                                                                          Areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 20% to 49% percent of total cover. These areas most commonly include single-family housing units.
    ## 5                                                                                                                                                                              Areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 50% to 79% of the total cover. These areas most commonly include single-family housing units.
    ## 6                                                                                                                                                          Highly developed areas where people reside or work in high numbers. Examples include apartment complexes, row houses and commercial/industrial. Impervious surfaces account for 80% to 100% of the total cover.
    ## 7                                                                                                                          Areas of bedrock, desert pavement, scarps, talus, slides, volcanic material, glacial debris, sand dunes, strip mines, gravel pits and other accumulations of earthen material. Generally, vegetation accounts for less than 15% of total cover.
    ## 8                                                                                                                                                             Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75% of the tree species shed foliage simultaneously in response to seasonal change.
    ## 9                                                                                                                                                  Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75% of the tree species maintain their leaves all year. Canopy is never without green foliage.
    ## 10                                                                                                                                                                        Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. Neither deciduous nor evergreen species are greater than 75% of total tree cover.
    ## 11                                                                                                                                        Alaska only areas dominated by shrubs less than 20 centimeters tall with shrub canopy typically greater than 20% of total vegetation. This type is often co-associated with grasses, sedges, herbs, and non-vascular vegetation.
    ## 12                                                                                                                        Areas dominated by shrubs; less than 5 meters tall with shrub canopy typically greater than 20% of total vegetation. This class includes true shrubs, young trees in an early successional stage or trees stunted from environmental conditions.
    ## 13                                                                                                                                                            Areas dominated by gramanoid or herbaceous vegetation, generally greater than 80% of total vegetation. These areas are not subject to intensive management such as tilling, but can be utilized for grazing.
    ## 14                                                                                                                                    Alaska only areas dominated by sedges and forbs, generally greater than 80% of total vegetation. This type can occur with significant other grasses or other grass like plants, and includes sedge tundra, and sedge tussock tundra.
    ## 15                                                                                                                                                                                                                                                             Alaska only areas dominated by fruticose or foliose lichens generally greater than 80% of total vegetation.
    ## 16                                                                                                                                                                                                                                                                                  Alaska only areas dominated by mosses, generally greater than 80% of total vegetation.
    ## 17                                                                                                                                     Areas of grasses, legumes, or grass-legume mixtures planted for livestock grazing or the production of seed or hay crops, typically on a perennial cycle. Pasture/hay vegetation accounts for greater than 20% of total vegetation.
    ## 18                                                                      Areas used for the production of annual crops, such as corn, soybeans, vegetables, tobacco, and cotton, and also perennial woody crops such as orchards and vineyards. Crop vegetation accounts for greater than 20% of total vegetation. This class also includes all land being actively tilled.
    ## 19                                                                                                                                                                                            Areas where forest or shrubland vegetation accounts for greater than 20% of vegetative cover and the soil or substrate is periodically saturated with or covered with water.
    ## 20                                                                                                                                                                                           Areas where perennial herbaceous vegetation accounts for greater than 80% of vegetative cover and the soil or substrate is periodically saturated with or covered with water.

``` r
levels(land_use)<-rat
```

Note that incoming NLCD data is the world web mercator projection
`"+init=epsg:3857"` this caused no end of grief in this weeks lab.
Therefore, we will transform the fisherSP2 object to be in this same
epsg projection to match the NLCD data using a comnination of the
commands `st_transform()` and `st_crs()`.

``` r
fisher_sf3 <-st_transform(fisher_sf2, st_crs(land_use)) 

tmap_mode("plot")
map <- tm_shape(land_use) +tm_raster(legend.show = FALSE)
map + tm_shape(fisher_sf3) +tm_dots()
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
plot(land_use)
plot(fisher_sf3, add=TRUE, type="p",cex=0.5)
```

![](README_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

Try adding a raster to our earlier mapview using + land_use. Note you
get a bunch of error messages, and it takes a LONNG time. But its a cool
way to explore raster data in mapview.

    mapview(fisherSP3, zcol="individual.local.identifier", legend = TRUE, cex=5, lwd=2, map.type = c("OpenStreetMap.DE", "Esri.WorldShadedRelief")) + land_use

## Exploring Landcover for RickyT

*adapted from Singer et al.*

We hypothesized that Fishers prefers forested wetlands over other
landuse classes, based on previous work by Scott LaPoint. Armed with our
NLCD landuse raster we create a layer called wet that is 1 for forested
wetlands (category 90) and 0 otherwise using the terra packages.

``` r
wet <- land_use == 90
names(wet) <- "wet"
#ext(wet) <- c(xmin=-8230000, xmax=-8210000, ymin=5270000, ymax=5280000)


##Lets zoom into Ricky T

rickyT.raster <- rast()
ext(rickyT.raster) <- c(xmin=1776730, xmax=1782940, ymin=2411165, ymax=2413945)
crs(rickyT.raster) <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"



plot(wet, ext = rickyT.raster)
plot(fisher_sf3,add=TRUE,  type="p", pch=12, cex = 0.5, ext = rickyT.raster)
```

    ## Warning in plot.sf(fisher_sf3, add = TRUE, type = "p", pch = 12, cex = 0.5, :
    ## ignoring all but the first attribute

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Wow, Ricky surely seems to preferentially select ‘wet’ forests. For our
first SSF we will focus just on this one covariate.

## Exploratory Analyses of Step Lengths and Turning Angles

Before fitting our SSF, we will conduct some explortatory analyses of
whether step lengtha and turning angles differ between wet and all other
landcover models. Our first step is to extract covariates for the start
point of each step using the amt::extract_covariates() function.

*from Singer et al*

Note that the function amt::extract \_ covariates takes an argument
where that indicates whether covariate values should be extracted at the
beginning or the end of a step (“both” can be used to extract the
covariate at the start and the end of a step). Depending on the target
process under investigation (habitat selection or movement), covariates
might be extracted at the end of the step (habitat selection process) or
at the start of the step (movement process).

If covariates are extracted at the end of the step, they are typically
included in the model as main effects, to answer questions of the type:
How do covariates influence where the animal moves?

In contrary, if covariates are extracted at the beginning of the step,
they are typically included in the model as an interaction with movement
characteristics (step length, log of the step length, or the cosine of
the turn angle), to test hypotheses of the type: Do animals move
faster/more directed, if they start in a given habitat?

Finally, covariate values at the start and the end of a step can also be
included in the model as an interaction with each other, to test
hypotheses of the type: Are animals more likely to stay in a given
habitat, if they are already in that habitat?

Here, we will….

``` r
#wet_layer <- as(wet, "Raster")

dat1 <- amt::make_track(dat_1, x, y, t, crs = 4326) %>%
     amt::transform_coords(st_crs(wet))

stps <- amt::track_resample(dat1, rate = minutes(10), tolerance = minutes(1)) %>%
  amt::filter_min_n_burst(min_n = 3) %>% amt::steps_by_burst() %>%
  amt::time_of_day(include.crepuscule = FALSE)

eda1 <- stps %>% amt::extract_covariates(wet, where ="start") %>% 
  mutate(wet = if_else(wet == "TRUE", 1, 0)) %>% 
  mutate(landuse = factor(wet, levels = c(0, 1), labels = c("other", "forested wetland")))
head(eda1)
```

    ## # A tibble: 6 × 14
    ##   burst_      x1_     x2_    y1_    y2_    sl_ direc…¹   ta_ t1_                
    ## *  <dbl>    <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dttm>             
    ## 1     16 1780410.  1.78e6 2.41e6 2.41e6 199.     -1.27 NA    2010-02-10 03:00:35
    ## 2     16 1780469.  1.78e6 2.41e6 2.41e6   6.01    1.70  2.96 2010-02-10 03:10:39
    ## 3     16 1780468.  1.78e6 2.41e6 2.41e6  20.8    -1.00 -2.70 2010-02-10 03:20:10
    ## 4     16 1780480.  1.78e6 2.41e6 2.41e6  11.2     2.72 -2.56 2010-02-10 03:30:24
    ## 5     18 1780471.  1.78e6 2.41e6 2.41e6   6.85   -1.87 NA    2010-02-10 05:00:43
    ## 6     18 1780469.  1.78e6 2.41e6 2.41e6   2.87    1.77 -2.64 2010-02-10 05:10:14
    ## # … with 5 more variables: t2_ <dttm>, dt_ <drtn>, tod_end_ <fct>, wet <dbl>,
    ## #   landuse <fct>, and abbreviated variable name ¹​direction_p

Next, we make some summary plots of step length, turning angle as a
function of day/night, and wet and other landcover types. We bundle them
together using the `cowplot` package.

``` r
p1 <- eda1 %>% dplyr::select(landuse, tod = tod_end_, sl_, ta_) %>%
  gather(key, val, -landuse, -tod) %>%
  filter(key == "sl_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ landuse, nrow = 2) +
  xlab("Step length [m]") + theme_light() +
  ylab("Density") +
  theme(legend.title = element_blank())
```

    ## Warning: attributes are not identical across measure variables;
    ## they will be dropped

``` r
p2 <- eda1 %>% dplyr::select(landuse, tod = tod_end_, sl_, ta_) %>%
  gather(key, val, -landuse, -tod) %>%
  filter(key == "ta_") %>%
  ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
  facet_wrap(~ landuse, nrow = 2) +
  xlab("Turn angle") + theme_light() +
  theme(legend.title = element_blank(),
  axis.title.y = element_blank())
```

    ## Warning: attributes are not identical across measure variables;
    ## they will be dropped

``` r
pg1 <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), rel_widths = c(1, 1))
```

    ## Warning: Removed 307 rows containing non-finite values (`stat_density()`).

``` r
leg <- get_legend(p1)
plot_grid(pg1, leg, rel_widths = c(1, 0.1))
```

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- --> If you want
to save the figure, use this

    #ggsave("Figures/fig_eda_1_animal.pdf", width = 20, height = 18, units = "cm")

This figure represents Exploratory data analysis of one individual
fisher, Ricky T (id: 1016): empirical distributions of step lengths
(first column) and turning angles (second column) are shown for forested
wetland (second row) and other habitats (first row) and for day and
night (colors).

We note that the underlying gamma distributions for the step lengths
vary by time of day, and, by landcover types such that it seems the
Fisher moves further at night and especially in forested wetlands.
Similarly, we see differences in the turning angles where there is much
stronger directional persistence during the night in forested wetlands.
We will see the coefficients for these differences in the fitting of the
SSF model.

## Fitting a Step Selection Function

*From Singer et al.*

To fit SSFs, the observed covariates associated with observed steps are
compared to covariates associated with random (or control) steps. Random
steps can be generated by either: (a) sampling from the observed turn
step‐length and turn angle distribution (resulting in a traditional SSF
al la Fortin et al. 2005), or (b) by fitting a parametric distribution
to the observed step lengths (either a negative‐exponential, a
half-normal, a log‐normal, or a gamma; see Avgar et al., 2016, Appendix
2) and turn angles (a von Mises; Duchesne et al., 2015).

As mentioned above, choosing b), an iSSF, is arguably less biased and
also provides the user with a mechanistic movement model that can be
used to simulate space use, and hence utilization distributions (Avgar
et al., 2016; Signer et al., 2017). Currently, amt only implements the
iSSFs with gamma and von Mises distributions.

### amt::random_steps

Before we fit the ssf, we must first generate random steps using
amt::random_steps(n=9) which choses 9 random points per Fisher location.
We then also extract covariates at the end points of each random step
(the wet covariate here), and add time of day and the log of step length
to each observation. We proceed by fitting a gamma distribution to the
step lengths and a von Mises distribution to the turn angles using
maximum likelihood (Agostinelli & Lund, 2017; Delignette‐Muller &
Dutang, 2015), and use these distributions to generate and pair nine
random steps with each observed step. The number of random steps effects
the estimation error; the more the steps, the lower the error, but the
higher the computational burden (Avgar et al., 2016).

``` r
m1 <-stps %>% amt::random_steps(n = 9) %>%
  amt::extract_covariates(wet) %>%
  mutate(wet = if_else(wet == "TRUE", 1, 0)) %>% 
  amt::time_of_day(include.crepuscule = FALSE) %>%
  mutate(unique_step = paste(burst_,step_id_,sep="_")) %>% 
  mutate(log_sl_ = log(sl_)) -> d1
```

To see what the random_steps() function did, take a look at the first 18
rows.

``` r
head(m1, n=18)
```

    ## # A tibble: 18 × 16
    ##    burst_      x1_      x2_      y1_     y2_    sl_      ta_ t1_                
    ##  *  <dbl>    <dbl>    <dbl>    <dbl>   <dbl>  <dbl>    <dbl> <dttm>             
    ##  1     16 1780469. 1780468. 2412224.  2.41e6   6.01  2.96    2010-02-10 03:10:39
    ##  2     16 1780469. 1780291. 2412224.  2.41e6 243.   -1.13    2010-02-10 03:10:39
    ##  3     16 1780469. 1780478. 2412224.  2.41e6  23.1   0.0958  2010-02-10 03:10:39
    ##  4     16 1780469. 1780494. 2412224.  2.41e6  51.7   2.34    2010-02-10 03:10:39
    ##  5     16 1780469. 1780258. 2412224.  2.41e6 229.   -1.48    2010-02-10 03:10:39
    ##  6     16 1780469. 1780547. 2412224.  2.41e6 146.    2.28    2010-02-10 03:10:39
    ##  7     16 1780469. 1780447. 2412224.  2.41e6  29.7  -2.59    2010-02-10 03:10:39
    ##  8     16 1780469. 1780470. 2412224.  2.41e6   4.09  0.00401 2010-02-10 03:10:39
    ##  9     16 1780469. 1780499. 2412224.  2.41e6  90.5   0.0342  2010-02-10 03:10:39
    ## 10     16 1780469. 1780278. 2412224.  2.41e6 195.   -1.68    2010-02-10 03:10:39
    ## 11     16 1780468. 1780480. 2412230.  2.41e6  20.8  -2.70    2010-02-10 03:20:10
    ## 12     16 1780468. 1780939. 2412230.  2.41e6 571.   -1.09    2010-02-10 03:20:10
    ## 13     16 1780468. 1780253. 2412230.  2.41e6 216.    1.51    2010-02-10 03:20:10
    ## 14     16 1780468. 1780468. 2412230.  2.41e6   3.41  2.99    2010-02-10 03:20:10
    ## 15     16 1780468. 1780505. 2412230.  2.41e6  61.5  -2.63    2010-02-10 03:20:10
    ## 16     16 1780468. 1780510. 2412230.  2.41e6  42.5  -1.46    2010-02-10 03:20:10
    ## 17     16 1780468. 1780436. 2412230.  2.41e6  47.9   0.626   2010-02-10 03:20:10
    ## 18     16 1780468. 1780472. 2412230.  2.41e6   4.78 -1.06    2010-02-10 03:20:10
    ## # … with 8 more variables: t2_ <dttm>, dt_ <drtn>, tod_end_ <fct>, case_ <lgl>,
    ## #   step_id_ <dbl>, wet <dbl>, unique_step <chr>, log_sl_ <dbl>

``` r
#str(m1)
```

The target variable case\_ is one for observed steps and zero for random
(or control) steps. Each step is paired with several (here 9) control
steps that form together a stratum (indicated by strat(step_id\_) in the
model formula). The function amt::random_steps automatically creates a
new column, step_id\_ , that identifies different strata.

Note that it is challenging to visualize the SSF point generation
process, and at this point, we note that the x1\_ y1\_ and x2\_ and y2\_
for the ‘random’ cases. I have tried to visualize the random steps here,
but I think because of the short duration, they do not seem as dramatic.

``` r
m1$caseF <-as.factor(m1$case_)
ggplot(m1, aes(x2_ ,y2_, colour = caseF)) + geom_point(aes(size = caseF, colour = caseF))
```

    ## Warning: Using size for a discrete variable is not advised.

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

### Fitting SSF models using amt::fit_issf()

Next we fit 3 different statistical Step dplyr::selection Function
models using the amt::fit_issf to fit a conditional logistic regression
model to the resulting data including movement‐related covariates with
the function amt::fit \_ issf (a wrapper to survival::clogit; Therneau &
Grambsch, 2000).

We included three different combinations of different main effects, and
their interactions, in 3 different SSF models below. These included the
environmental covariate wet, the step length, and and the log of the
step length (log_sl\_ ) as modifiers of the shape parameter of the
underlying gamma distribution. The estimated coefficient of sl\_ and
log_sl\_ can be used to adjust the tentative shape estimate (i.e., the
estimate of the shape parameter using the observed step lengths) of the
underlaying gamma distribution for the step lengths. We also include
interactions between wet and tod\_ ,a factor with two levels—day (the
reference category) and night, and between tod\_ and log_sl\_ .These
interactions are included to the test the hypotheses that habitat
dplyr::selection and displacement rate, respectively, differ between day
and night.

We could have also included cosines of the turning angles and their
interaction with day. This choice would modify the concentration
parameter of the underlying von Mises distribution for the turning
angles and allow the degree of directional persistence to depend on time
of day; the data summarized in Figure 1 suggest that this could be a
sensible choice. For the sake of simplicity, however, we have assumed we
have correctly modeled the degree of directional persistence and that it
does not differ between day and night.

``` r
m3 <- d1 %>% amt::fit_issf(case_ ~ wet + sl_ + wet:tod_end_+ sl_:tod_end_ +strata(unique_step))
m2 <- d1 %>% amt::fit_issf(case_ ~ wet + log_sl_ + wet:tod_end_+ log_sl_:tod_end_ + strata(unique_step))
m1 <- d1 %>% amt::fit_issf(case_ ~ wet + log_sl_ + sl_ + wet:tod_end_+ log_sl_:tod_end_ + sl_:tod_end_ + strata(unique_step))
```

Model 3 represents the habitat hypothesis that Fisher resource
dplyr::selection differs for wet forests, and, differs between day and
night (represented by the interaction between wet:tod_end\_). Movement
hypotheses addressed as that step lenghts also differ between times of
day, for example. Here, the absolute step length, sl\_, is modeled.

Model 2 represents a similar hypothesis, but instead, with the log of
step length. Finally, Model 1 includes both the untransformed step
length and the log of step length.

### Model Selection

We can conduct model dplyr::selection using AIC on these 3 models.
First, I note however, that the Log-Likelihood for a conditional
logistic regression model is not directly comparable to that of a
traditional logistic regression model. So we cannot test, using AIC,
whether conditioning on availability at every time step ‘improves’ model
fit relative to a ‘naive’ logistic regression model. We will demonstrate
this below when we explore mixed-effects SSF models.

``` r
AIC(m1$model, m2$model, m3$model)
```

    ##          df      AIC
    ## m1$model  6 5408.965
    ## m2$model  4 5415.015
    ## m3$model  4 5405.042

So clearly, model1 is a bit better than model 2 and 3. Lets look at
Model 1, and summarize/extract the coefficients.

``` r
summary(m1)
```

    ## Call:
    ## coxph(formula = Surv(rep(1, 11870L), case_) ~ wet + log_sl_ + 
    ##     sl_ + wet:tod_end_ + log_sl_:tod_end_ + sl_:tod_end_ + strata(unique_step), 
    ##     data = data, method = "exact")
    ## 
    ##   n= 11870, number of events= 1187 
    ## 
    ##                            coef exp(coef)  se(coef)      z Pr(>|z|)    
    ## wet                    0.923750  2.518718  0.253276  3.647 0.000265 ***
    ## log_sl_               -0.024776  0.975528  0.100077 -0.248 0.804465    
    ## sl_                   -0.003418  0.996588  0.001603 -2.133 0.032959 *  
    ## wet:tod_end_night     -0.445828  0.640294  0.269524 -1.654 0.098100 .  
    ## log_sl_:tod_end_night  0.030286  1.030749  0.108769  0.278 0.780675    
    ## sl_:tod_end_night      0.004412  1.004421  0.001665  2.649 0.008074 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                       exp(coef) exp(-coef) lower .95 upper .95
    ## wet                      2.5187     0.3970    1.5332    4.1378
    ## log_sl_                  0.9755     1.0251    0.8018    1.1869
    ## sl_                      0.9966     1.0034    0.9935    0.9997
    ## wet:tod_end_night        0.6403     1.5618    0.3775    1.0859
    ## log_sl_:tod_end_night    1.0307     0.9702    0.8329    1.2757
    ## sl_:tod_end_night        1.0044     0.9956    1.0011    1.0077
    ## 
    ## Concordance= 0.556  (se = 0.011 )
    ## Likelihood ratio test= 69.37  on 6 df,   p=5e-13
    ## Wald test            = 64.28  on 6 df,   p=6e-12
    ## Score (logrank) test = 66.84  on 6 df,   p=2e-12

``` r
s <- summary(m1$model)$coefficients
s
```

    ##                               coef exp(coef)    se(coef)          z
    ## wet                    0.923750033 2.5187180 0.253276271  3.6472032
    ## log_sl_               -0.024776392 0.9755280 0.100077466 -0.2475721
    ## sl_                   -0.003418104 0.9965877 0.001602801 -2.1325815
    ## wet:tod_end_night     -0.445828366 0.6402937 0.269523662 -1.6541344
    ## log_sl_:tod_end_night  0.030285652 1.0307489 0.108769037  0.2784400
    ## sl_:tod_end_night      0.004411608 1.0044214 0.001665406  2.6489678
    ##                           Pr(>|z|)
    ## wet                   0.0002651102
    ## log_sl_               0.8044654720
    ## sl_                   0.0329590710
    ## wet:tod_end_night     0.0981002044
    ## log_sl_:tod_end_night 0.7806746112
    ## sl_:tod_end_night     0.0080738028

*from Singer et al. 2019* Inspecting the fitted model, we make the
following observations: (a) There is evidence to suggest that the animal
prefers forested wetlands over other landuse classes, (b) there is no
difference in RickyT’s preference for wetlands between day and night,
(c) there is evidence to modify the shape of the gamma distribution fit
to the observed step lengths (through the log of the step length), and
(d) the modification of the shape parameter should be done separately
for day and night, indicating that expected movement speeds differ
between day and night.

### Examining Movement Statistics

Besides inspecting the coefficients and their standard errors, we can
calculate movement statistics and examine the turning angle. We begin by
retrieving the tentative parameter estimates (i.e., the estimated
parameters before correcting for habitat dplyr::selection; see Avgar et
al. (2016) for more details) for the gamma distribution of the
step‐length distribution:

``` r
coef(m1)
```

    ##                   wet               log_sl_                   sl_ 
    ##           0.923750033          -0.024776392          -0.003418104 
    ##     wet:tod_end_night log_sl_:tod_end_night     sl_:tod_end_night 
    ##          -0.445828366           0.030285652           0.004411608

``` r
plot_sl(m1)
```

![](README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
sl_m1<-sl_distr(m1)
## note we will extract the scale and shape parameters of the gamma distribution for later simulating the UD
scale <- sl_m1$params$scale
shape <- sl_m1$params$shape
```

Note that these scale and shape parameters are for the intercet step
length distribution for the model, which is during the DAY. See the
coef(m1) output above for a reminder, and the figure we made above.

``` r
ggplot(d1, aes(tod_end_, sl_)) +geom_violin()
```

![](README_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

``` r
sl_data <- as_tibble(d1)
sl_data %>% group_by(tod_end_) %>% summarize(median = median(sl_))
```

    ## # A tibble: 2 × 2
    ##   tod_end_ median
    ##   <fct>     <dbl>
    ## 1 day        69.1
    ## 2 night      77.3

Note the interpretation here is that at night, fishers move \~ 105m per
step length, 7 meters further than during the day.

## Simulating the Utilization Distributions

In a final step, we simulated space‐use from the fitted model m1 to
obtain a model‐based estimate of the animal’s utilization distribution
(UD; Avgar et al., 2016; Signer et al., 2017). Generally, two types of
UDs can be simulated: the transient UD and the steady‐state UD. The
transient UD describes the expected space‐use distribution of the animal
within a short time period and is hence conditional on the starting
position. The steady‐state UD describes the expected space‐use
distribution of the animal in the long‐term. In order to simulate UDs,
one has to ensure that the animals stay within the study domain.

We see three possible methods for achieving this goal–all implemented in
amt: (a) use a covariate that attracts the animal toward one or more
centers of activity (e.g., the squared distance to the mean of all
coordinates), (b) use a very large landscape, or (c) use a wrapped
landscape (torus).

Here, we illustrate the simulation of steady‐state and transient UDs.
For the steady‐state UD, we simulate from the first observed location
107 time steps on a toroid landscape, once for day and once for night.
For the transient UD, we are interested in the UD up to 10 hr after last
observation, we therefore simulated 72 steps (at a 10 min sampling rate)
5 × 103 times.

First, we crop out wet for a smaller area around just Ricky T using the
amt::bbox() function which describes a bounding box.

``` r
wet_c <- crop(wet, amt::bbox(dat1, spatial = TRUE, buff = 1e3))
```

First, we recall that since we have different model coefficients for day
and night, we need to fit the movement and habitat kernels for day and
night.

Thus, we estimte the daytime movement kernel using the
amt::movement_kernel(). This calculates a movement kernel from a fitted
(i)SSF. The method is currently only implemented for the gamma
distribution. Note, we use the tentative scale estimate and the shape
estimate adjusted for day.

``` r
mk <- amt::movement_kernel(scale, shape, wet_c)
plot(mk)
```

Second, we estimate the habitat kernel which is calculated more or less
like a traditional RSF model by multiplying resources with their
corresponding coefficients from the fitted (i)SSF. That is for each
pixel we calculate the estimated dplyr::selection coefficients times the
resources and exponentiate the product.

``` r
#wet_c_rast <- as(wet_c, "Raster") #habitat_kernel needs the resource argument to be a raster layer or raster stack (wet_c is a formal class SpatRaster from terra)
hk <- amt::habitat_kernel(list(wet = coef(m1)["wet"]), wet_c)
plot(hk)
```

Next, we fit the simulated steady-state UD based on the movement model
embedded within the SSF using amt:: simulate_ud(), which simulates a
utilization distribution (UD) from a fitted Step-dplyr::selection
Function. We time it because we are into that sort of thing.

``` r
system.time(ssud_day <- amt::simulate_ud(
  mk, hk,
  as.numeric(stps[1, c("x1_", "y1_")]),
   n = 1e7))
 plot(ssud_day)
```

Finally, we simulate the transient UD which is based on some starting
location and duration (72 steps here) using the amt:simulate_tud(). This
is a conviencience wrapper arround simulate_ud to simulate transition
UDs (i.e., starting at the same position many times and only simulate
for a short time). In order to simulate the transient UD we have to
repeatedly simulate short tracks starting at the same point, and then
sum individual UDs and normalize, which we do with the function
amt::simulate \_ tud.

``` r
system.time(tud_day <- amt::simulate_tud(mk, hk, as.numeric(stps[150, c("x1_", "y1_")]), n = 72, n_rep = 5e3))
plot(tud_day)
```

And finally, we bundle them together in a nice figure

``` r
pllog <- list(
  geom_raster(),
   coord_equal(),
  scale_fill_continuous(low = "white", high = "red", tran = "log10", na.value = "white"),
   scale_y_continuous(expand = c(0, 0)),
   scale_x_continuous(expand = c(0, 0)),
   theme_light(),
   theme(legend.position = "none"))

 pl <- list(
   geom_raster(),
   coord_equal(),
   scale_fill_continuous(low = "white", high = "red", na.value = "white"),
 scale_y_continuous(expand = c(0, 0)),
   scale_x_continuous(expand = c(0, 0)),
   theme_light(),
   theme(legend.position = "none"))

r1 <- as.data.frame(st_as_stars(mk)) #converting raster mk to a stars object to get x and y coordinates
p1 <- ggplot(r1, aes(x, y, fill = d)) + pllog + ggtitle("Movement kernel (day)")


r2 <- as.data.frame(st_as_stars(hk))
p2 <- ggplot(r2, aes(x, y, fill = layer)) + pl + ggtitle("Habitat kernel (day)")

r1 <- as.data.frame(st_as_stars(tud_day))
 p3 <- ggplot(r1, aes(x, y, fill = layer)) + pllog + ggtitle("Transient UD (day)")

r1 <- as.data.frame(st_as_stars(ssud_day))
p5 <- ggplot(r1, aes(x, y, fill = layer)) + pl + ggtitle("Steady state UD (day)")

cowplot::plot_grid(p1, p2, p3, p5, ncol = 2, labels = "AUTO")
```

This figure shows Simulated utilization distributions for 1 fisher,
RickyT, during the DAY. To obtain simulated Utilization Distributions
(UD), a movement kernel (panel a) and a habitat kernel (panel b) are
needed. The movement kernel is always placed at the current position of
the animal. The next step of the animal is then sampled with probability
proportional to the product of two kernels. Expected differences in
movement speeds between night and day are reflected in the transient UD
(panels c and e) and to a lesser extend in steady‐state UD (panels d and
f). Note, for better visualization, fills were log10 transformed for
panels a, c, and e

Note that the transient UD’s are slightly different from the paper
because I started in a slightly different starting location.

Finally, if we wanted to save them all:

    ggsave("fig_one_animal1.pdf", height = 20, width = 24, units = "cm")

Also note that the paper in Ecology and Evolution shows calculation of
the UD’s during the night, which I have not been able to reproduce here
because the `amt` package no longer has the adjust_shape() function.
Welcome to the cutting edge if you want to learn how to adjust movement
parameters. I can think of dividing the data into two subsets, day and
night.

# Fitting SSF Models to Multiple Animals

We will conduct two-stage modeling of SSF’s fit to individual Fishers in
the same study area.

We start again with the same data set (dat), containing data from six
individual fishers. This time we are interested in quantifying
among‐animal variability in the dplyr::selection coefficients. We
proceed using nearly all the same steps as in the first example, but
with a different data structure: data\_ frames with list columns (Müller
& Wickham, 2018). List columns are best thought of as regular columns of
a data\_ frame that are R lists and can contain any objects (in our case
tracks and fitted models). The purrr::nest command can be used to nest
data into a list column (Henry & Wickham, 2017).

``` r
dat <- read_csv("Data/Martes pennanti LaPoint New York.csv") %>%
 filter(!is.na(`location-lat`)) %>%
  dplyr::select(x = `location-long`, y = `location-lat`,
              t = `timestamp`, id = `tag-local-identifier`) %>%
   filter(id %in% c(1465, 1466, 1072, 1078, 1016, 1469))

 dat_all <- dat %>% nest(-id)
dat_all$sex <- c("f", "f", "f", "m", "m", "m")
dat_all
```

    ## # A tibble: 6 × 3
    ##      id data                 sex  
    ##   <dbl> <list>               <chr>
    ## 1  1072 <tibble [1,349 × 3]> f    
    ## 2  1465 <tibble [3,004 × 3]> f    
    ## 3  1466 <tibble [1,501 × 3]> f    
    ## 4  1078 <tibble [1,638 × 3]> m    
    ## 5  1469 <tibble [2,436 × 3]> m    
    ## 6  1016 <tibble [8,958 × 3]> m

dat_all is now a nested data frame with 6 rows (one for each individual)
and two columns. In the first column the animal id is given, and in the
second column (by default named data) the relocations of the
corresponding animal are saved.

We can now apply the steps as before for all animals. We first create a
track for each animal and transform the coordinate reference system
using the function amt::transform \_ coords.

``` r
 dat_all <- dat_all %>%
  mutate(trk = map(data, function(d) {
     amt::make_track(d, x, y, t, crs = 4326) %>%
       amt::transform_coords(st_crs(land_use))}))
```

And summarize sampling rate

``` r
dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% unnest(cols = c(sr))
```

    ## # A tibble: 6 × 10
    ##      id   min    q1 median  mean    q3   max    sd     n unit 
    ##   <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <int> <chr>
    ## 1  1072 3.92   9.80  10.0  20.7  10.3  1650.  95.2  1348 min  
    ## 2  1465 0.400  1.97   2.03  8.93  3.98 1080.  48.1  3003 min  
    ## 3  1466 0.417  1.97   2.07 15.7   4.08 2167. 104.   1500 min  
    ## 4  1078 1.35   9.78  10.0  21.7  10.3  2111.  87.2  1637 min  
    ## 5  1469 0.417  1.97   2.17 13.3   5.42 2889.  90.2  2435 min  
    ## 6  1016 0.100  1.93   2.03  8.04  2.57 1209.  44.0  8957 min

This time we see that some individuals have a 2 min sample rate and
others a 10 min one. Thus, we decided to resample the tracks to the same
sampling rate of 10 min (noting that (i)SSF inference is scale
dependent; Signer et al., 2017) using amt::track_resample.

## NLCD Extraction

Next, we are going to test for SSF as a function of 4 different
landcover types now, more than just the wet forests. So we will create a
reclassification matrix that lumps everything into 5 categories, with
water/wet forests as the reference catetgory: water and wetland forests,
developed open spaces, other developed areas, forests and shrubs, and
crops.

First, we remind ourselves of the Raster Attribute Table (RAT) we
created above.

``` r
rat
```

    ##    ID                        Class   Color
    ## 1  11                   Open Water #5475A8
    ## 2  12           Perennial Ice/Snow #FFFFFF
    ## 3  21        Developed, Open Space #E8D1D1
    ## 4  22     Developed, Low Intensity #E29E8C
    ## 5  23  Developed, Medium Intensity #ff0000
    ## 6  24     Developed High Intensity #B50000
    ## 7  31 Barren Land (Rock/Sand/Clay) #D2CDC0
    ## 8  41             Deciduous Forest #85C77E
    ## 9  42             Evergreen Forest #38814E
    ## 10 43                 Mixed Forest #D4E7B0
    ## 11 51                  Dwarf Scrub #AF963C
    ## 12 52                  Shrub/Scrub #DCCA8F
    ## 13 71         Grassland/Herbaceous #FDE9AA
    ## 14 72             Sedge/Herbaceous #D1D182
    ## 15 73                      Lichens #A3CC51
    ## 16 74                         Moss #82BA9E
    ## 17 81                  Pasture/Hay #FBF65D
    ## 18 82             Cultivated Crops #CA9146
    ## 19 90               Woody Wetlands #C8E6F8
    ## 20 95 Emergent Herbaceous Wetlands #64B3D5
    ##                                                                                                                                                                                                                                                                                                                                                                Description
    ## 1                                                                                                                                                                                                                                                                                           Areas of open water, generally with less than 25% cover of vegetation or soil.
    ## 2                                                                                                                                                                                                                                                                  Areas characterized by a perennial cover of ice and/or snow, generally greater than 25% of total cover.
    ## 3  Areas with a mixture of some constructed materials, but mostly vegetation in the form of lawn grasses. Impervious surfaces account for less than 20% of total cover. These areas most commonly include large-lot single-family housing units, parks, golf courses, and vegetation planted in developed settings for recreation, erosion control, or aesthetic purposes.
    ## 4                                                                                                                                                                          Areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 20% to 49% percent of total cover. These areas most commonly include single-family housing units.
    ## 5                                                                                                                                                                              Areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 50% to 79% of the total cover. These areas most commonly include single-family housing units.
    ## 6                                                                                                                                                          Highly developed areas where people reside or work in high numbers. Examples include apartment complexes, row houses and commercial/industrial. Impervious surfaces account for 80% to 100% of the total cover.
    ## 7                                                                                                                          Areas of bedrock, desert pavement, scarps, talus, slides, volcanic material, glacial debris, sand dunes, strip mines, gravel pits and other accumulations of earthen material. Generally, vegetation accounts for less than 15% of total cover.
    ## 8                                                                                                                                                             Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75% of the tree species shed foliage simultaneously in response to seasonal change.
    ## 9                                                                                                                                                  Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75% of the tree species maintain their leaves all year. Canopy is never without green foliage.
    ## 10                                                                                                                                                                        Areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. Neither deciduous nor evergreen species are greater than 75% of total tree cover.
    ## 11                                                                                                                                        Alaska only areas dominated by shrubs less than 20 centimeters tall with shrub canopy typically greater than 20% of total vegetation. This type is often co-associated with grasses, sedges, herbs, and non-vascular vegetation.
    ## 12                                                                                                                        Areas dominated by shrubs; less than 5 meters tall with shrub canopy typically greater than 20% of total vegetation. This class includes true shrubs, young trees in an early successional stage or trees stunted from environmental conditions.
    ## 13                                                                                                                                                            Areas dominated by gramanoid or herbaceous vegetation, generally greater than 80% of total vegetation. These areas are not subject to intensive management such as tilling, but can be utilized for grazing.
    ## 14                                                                                                                                    Alaska only areas dominated by sedges and forbs, generally greater than 80% of total vegetation. This type can occur with significant other grasses or other grass like plants, and includes sedge tundra, and sedge tussock tundra.
    ## 15                                                                                                                                                                                                                                                             Alaska only areas dominated by fruticose or foliose lichens generally greater than 80% of total vegetation.
    ## 16                                                                                                                                                                                                                                                                                  Alaska only areas dominated by mosses, generally greater than 80% of total vegetation.
    ## 17                                                                                                                                     Areas of grasses, legumes, or grass-legume mixtures planted for livestock grazing or the production of seed or hay crops, typically on a perennial cycle. Pasture/hay vegetation accounts for greater than 20% of total vegetation.
    ## 18                                                                      Areas used for the production of annual crops, such as corn, soybeans, vegetables, tobacco, and cotton, and also perennial woody crops such as orchards and vineyards. Crop vegetation accounts for greater than 20% of total vegetation. This class also includes all land being actively tilled.
    ## 19                                                                                                                                                                                            Areas where forest or shrubland vegetation accounts for greater than 20% of vegetative cover and the soil or substrate is periodically saturated with or covered with water.
    ## 20                                                                                                                                                                                           Areas where perennial herbaceous vegetation accounts for greater than 80% of vegetative cover and the soil or substrate is periodically saturated with or covered with water.

``` r
rcl <- cbind(c(11, 12, 21:24, 31, 41:43, 51:52, 71:74, 81:82, 90, 95),c(1, 1, 2, 3, 3, 3, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 8, 8, 1, 1))
# water, dev open, dev, barren, forest, shrub and herb, crops, wetlands
# 1: water, wetlands
# 2: developed (open)
# 3: developed (other)
# 5: forest, herbaceous
 # 8: crops
lu <- classify(land_use, rcl, right = NA)
names(lu) <- "landuse"
plot(lu)
plot(fisher_sf3, add=TRUE)
```

    ## Warning in plot.sf(fisher_sf3, add = TRUE): ignoring all but the first attribute

![](README_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Next, we resample the tracks from 6 individual fishers to 10 minute
steps, filter them into consecutive bursts of 3 or more 10 minute
locations, create random steps (which by default goes back to the last
setting used, which was n=9), extract covariate values from both the
start and the end, and convert the landuse covariate to a factor all in
one step.

``` r
#lu2 <- as(lu, "Raster") #changing SpatRaster object lu into a Raster Layer for use in the command amt::extract_covariates which requires a raster layer

 m1 <- dat_all %>%
   mutate(steps = map(trk, function(x) {
     x %>% amt::track_resample(rate = minutes(10), tolerance = seconds(120)) %>%
      amt::filter_min_n_burst() %>%
       amt::steps_by_burst() %>% amt::random_steps() %>%
       amt::extract_covariates(lu, where = "both") %>% 
       mutate(unique_step = paste(burst_,step_id_,sep="_")) %>% 
       mutate(landuse_end = factor(landuse_end))
     }))
```

Note we got some NaNs produced errors - 12 of them to be precise. We’d
have to dig into exactly why, but I would not be TOO concerned in such a
big dataset. But - my guess is its because some available locations
ended off the map.

We then fit a simple SSF model to each individual. The main difference
to the previous example here, is that the all the steps from above are
wrapped into one mutate call. This call creates a new column to dat \_
all called ssf. This is a list column and each entry in this column
contains a fitted SSF.

``` r
m4 <- m1 %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ landuse_end +
                                                           strata(unique_step))))
m4
```

    ## # A tibble: 6 × 6
    ##      id data                 sex   trk                    steps      fit       
    ##   <dbl> <list>               <chr> <list>                 <list>     <list>    
    ## 1  1072 <tibble [1,349 × 3]> f     <trck_xyt [1,349 × 3]> <rndm_stp> <fit_clgt>
    ## 2  1465 <tibble [3,004 × 3]> f     <trck_xyt [3,004 × 3]> <rndm_stp> <fit_clgt>
    ## 3  1466 <tibble [1,501 × 3]> f     <trck_xyt [1,501 × 3]> <rndm_stp> <fit_clgt>
    ## 4  1078 <tibble [1,638 × 3]> m     <trck_xyt [1,638 × 3]> <rndm_stp> <fit_clgt>
    ## 5  1469 <tibble [2,436 × 3]> m     <trck_xyt [2,436 × 3]> <rndm_stp> <fit_clgt>
    ## 6  1016 <tibble [8,958 × 3]> m     <trck_xyt [8,958 × 3]> <rndm_stp> <fit_clgt>

m4 is still a data frame with one new column: ssf that is again a list
column with a fitted SSF. Lets look at 1 of these models for individual
1

``` r
m4$fit[[1]]$model
```

    ## Call:
    ## survival::clogit(formula, data = data, ...)
    ## 
    ##                  coef exp(coef) se(coef)      z        p
    ## landuse_end2 -0.58218   0.55868  0.19030 -3.059 0.002219
    ## landuse_end3 -2.31557   0.09871  0.47367 -4.889 1.02e-06
    ## landuse_end5  0.52918   1.69754  0.14008  3.778 0.000158
    ## landuse_end8 -1.08639   0.33743  0.41658 -2.608 0.009111
    ## 
    ## Likelihood ratio test=173  on 4 df, p=< 2.2e-16
    ## n= 13112, number of events= 1192

Lets look at the marginal, population level averaged coefficient for
landuse across our individual animals.

``` r
d2 <- m4 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
   dplyr::select(id, sex, coef) %>% unnest(cols = coef) %>%
  mutate(id = factor(id)) %>% group_by(term) %>%
  summarize(
     mean = mean(estimate),
    ymin = mean - 1.96 * sd(estimate),
     ymax = mean + 1.96 * sd(estimate))

d2$x <- 1:nrow(d2)
d2
```

    ## # A tibble: 4 × 5
    ##   term           mean   ymin   ymax     x
    ##   <chr>         <dbl>  <dbl>  <dbl> <int>
    ## 1 landuse_end2 -0.945 -2.06   0.169     1
    ## 2 landuse_end3 -2.06  -2.98  -1.15      2
    ## 3 landuse_end5  0.113 -0.625  0.850     3
    ## 4 landuse_end8 -1.32  -2.49  -0.156     4

Finally, we make some fancy figures

``` r
p1 <- m4 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model, conf.int = T))) %>%
  dplyr::select(id, sex, coef) %>% 
  unnest(cols = coef) %>% 
  mutate(id = factor(id)) %>%
  ggplot(., aes(x = term, y = estimate, group = id, col = id, pch = sex)) +
  geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin, ymax = ymax), data = d2, inherit.aes = FALSE,fill = "grey90") +
  geom_segment(mapping = aes(x = x - .4, xend = x + .4,y = mean, yend = mean), data = d2, inherit.aes = FALSE, size = 1) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.7), size = 0.8) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Habitat", y = "Relative selection Strength") + theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops"))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.

``` r
p1
```

![](README_files/figure-gfm/unnamed-chunk-39-1.png)<!-- --> From here,
it is easy to investigate coefficients for several animals and look at
population‐level effects. The results suggest that there are some
general population‐level trends. All fishers seem to prefer wetland
forests and natural areas relative to developed areas (of either type),
whereas considerable among‐animal variability in the coefficients for
crops makes it difficult to draw firm conclusions about this landuse
type. Lastly, there seems to be little differentiation based on sex.

And again, if we want to save we use:

    ggsave("img/fig_all_animals.pdf", width = 24, height = 12, units = "cm")

# Mixed-effect cLogit Models

In the next set of exercises, I will build on the SSF case study of
Singer et al., and conduct mixed-effects clogit analyses accounting for
individual fishers as random effects similar to Lab 7. Conceptually,
adding random effects to conditional logistic regression models was
challenging because there is no intercept. Here are the first 2 papers
that figured out how to add a random intercept for each individual
animal (e.g.), however, it did so in MATLAB. So, its mostly inaccessible
to biologists.

Craiu, R. V., T. Duchesne, D. Fortin, and S. Baillargeon. 2011.
Conditional Logistic Regression With Longitudinal Follow-up and
Individual-Level Random Coefficients: A Stable and Efficient Two-Step
Estimation Method. Journal of Computational and Graphical Statistics
20:767-784.

Duchesne, T., D. Fortin, and N. Courbin. 2010. Mixed conditional
logistic regression for habitat dplyr::selection studies. Journal of
Animal Ecology 79:548-555.

Since these initial papers, however, there have been a few big
breakthrough’s lately with the mclogit package
<http://cran.r-project.org/web/packages/mclogit/mclogit.pdf> or the
coxme package here
<http://cran.r-project.org/web/packages/coxme/coxme.pdf> I just played
around with both of these packages and they are actually.

We will first need to ‘unpack’ the nested data frame from section 3
above into an expanded dataframe using the unnest command. Then we will
progress through a set of 3-4 different models and compare model
interpretations to naive GLM models of the same kind.

``` r
fisher6 <- dat_all %>%
  mutate(steps = map(trk, function(x) {
    x %>% amt::track_resample(rate = minutes(10), tolerance = seconds(120)) %>%
      amt::filter_min_n_burst() %>%
      amt::steps_by_burst() %>% amt::random_steps() %>%
      amt::extract_covariates(lu, where = "both") %>%
      mutate(unique_step = paste(burst_,step_id_,sep="_")) %>% 
      mutate(landuse_end = factor(landuse_end))
  })) %>%
  dplyr::select(id, steps) %>%
  unnest()
```

    ## Warning: `cols` is now required when using unnest().
    ## Please use `cols = c(steps)`

``` r
fisher6
```

    ## # A tibble: 58,377 × 16
    ##       id burst_      x1_    x2_    y1_    y2_    sl_     ta_ t1_                
    ##    <dbl>  <dbl>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dttm>             
    ##  1  1072      4 1784400. 1.78e6 2.41e6 2.41e6  83.3   2.79   2011-02-11 19:50:58
    ##  2  1072      4 1784400. 1.78e6 2.41e6 2.41e6  34.8  -0.619  2011-02-11 19:50:58
    ##  3  1072      4 1784400. 1.78e6 2.41e6 2.41e6   1.37 -0.942  2011-02-11 19:50:58
    ##  4  1072      4 1784400. 1.78e6 2.41e6 2.41e6  75.3  -1.78   2011-02-11 19:50:58
    ##  5  1072      4 1784400. 1.78e6 2.41e6 2.41e6 201.   -2.85   2011-02-11 19:50:58
    ##  6  1072      4 1784400. 1.78e6 2.41e6 2.41e6 148.   -2.24   2011-02-11 19:50:58
    ##  7  1072      4 1784400. 1.78e6 2.41e6 2.41e6   6.04  2.32   2011-02-11 19:50:58
    ##  8  1072      4 1784400. 1.78e6 2.41e6 2.41e6  60.4  -0.0155 2011-02-11 19:50:58
    ##  9  1072      4 1784400. 1.78e6 2.41e6 2.41e6   1.06 -0.448  2011-02-11 19:50:58
    ## 10  1072      4 1784400. 1.78e6 2.41e6 2.41e6 113.    1.82   2011-02-11 19:50:58
    ## # … with 58,367 more rows, and 7 more variables: t2_ <dttm>, dt_ <drtn>,
    ## #   case_ <lgl>, step_id_ <dbl>, landuse_start <dbl>, landuse_end <fct>,
    ## #   unique_step <chr>

Next, we will add a String variable for landuse_end recalling that,
earlier, we defined the NLCD landcover according to: 1: water, wetlands
2: developed (open) 3: developed (other) 5: forest, herbaceouse 8: crops

``` r
head(fisher6$landuse_end)
```

    ## [1] 5 5 5 5 5 5
    ## Levels: 1 2 3 5 8

``` r
fisher6$landuseName = ifelse(fisher6$landuse_end == 1, "Wet Forests", 
                             ifelse(fisher6$landuse_end == 2, "Developed Open", 
                             ifelse(fisher6$landuse_end == 3, "Developed Other", 
                             ifelse(fisher6$landuse_end == 5, "Natural", "Crops"))))
table(fisher6$landuseName, fisher6$landuse_end)
```

    ##                  
    ##                       1     2     3     5     8
    ##   Crops               0     0     0     0  1858
    ##   Developed Open      0  5411     0     0     0
    ##   Developed Other     0     0  4877     0     0
    ##   Natural             0     0     0 29972     0
    ##   Wet Forests     16259     0     0     0     0

## Fit a naive GLM

First, we will fit a ‘naive’ GLM only focusing on the habitat processes,
that is, habitat dplyr::selection for the landcover covaraites and
compare them to the coefficients from the SSF fit to each individual
Fisher and their two-stage population-level averages.

``` r
model1 <- glm(case_~ I(landuse_end), data=fisher6,family=binomial(link="logit"))
## I commented out these next few versions of the models fit to landuseName to make comparisons to the previously fit 6 fisher two-step models more comparable, though we have to then keep track of which landovers 2, 3, 5, and 8 are. 
#model1 <- glm(case_~ I(landuseName), data=fisher6,family=binomial(link="logit"))
#model1 <- glm(case_~ I(landuseName=="Developed Open") + I(landuseName=="Developed Other") +I(landuseName=="Natural")+I(landuseName=="Crops"), data=fisher6,family=binomial(link="logit"))
summary(model1)
```

    ## 
    ## Call:
    ## glm(formula = case_ ~ I(landuse_end), family = binomial(link = "logit"), 
    ##     data = fisher6)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -0.4700  -0.4700  -0.4577  -0.3727   2.6982  
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -2.20332    0.02621 -84.079  < 2e-16 ***
    ## I(landuse_end)2 -0.42902    0.06033  -7.111 1.15e-12 ***
    ## I(landuse_end)3 -1.41034    0.09333 -15.112  < 2e-16 ***
    ## I(landuse_end)5  0.05581    0.03230   1.728    0.084 .  
    ## I(landuse_end)8 -1.11394    0.12895  -8.638  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 35567  on 58376  degrees of freedom
    ## Residual deviance: 35024  on 58372  degrees of freedom
    ## AIC: 35034
    ## 
    ## Number of Fisher Scoring iterations: 6

This model gives us an intercept, which is interpreted as wet-forests,
and we note that we also did not get an intercept in the SSF model
above. Lets now compare the coefficients to above:

``` r
coef(model1)
```

    ##     (Intercept) I(landuse_end)2 I(landuse_end)3 I(landuse_end)5 I(landuse_end)8 
    ##     -2.20332153     -0.42902304     -1.41033755      0.05581288     -1.11393667

``` r
naive_glm <- broom::tidy(model1) %>% 
  filter(!term=="(Intercept)") 

figNaive <- naive_glm %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate - 1.96*std.error, ymax = estimate +1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
    theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))
figNaive
```

![](README_files/figure-gfm/unnamed-chunk-43-1.png)<!-- --> First, we
recall we are comparing the two-stage averaged SSF coefficients from the
fisher SSF model above extracted in the object d2. This is simply the
arithmetic mean of the 6 individual coefficients. Second, note that the
ORDER of landcover categories on the X axis are changed now, which is
slightly annoying, and tough to remedy. But keep that in mind when
looking at the next figure.

``` r
fig5 <- plot_grid(p1, figNaive)
fig5
```

![](README_files/figure-gfm/unnamed-chunk-44-1.png)<!-- -->

We note that there are some similarities between coefficients, but,
differences.

Developed Open is -0.49 from the Naive Logit, and -1.00 from the SSF
Developed Other is -1.82 from the Naive Logit, and -2.20 from the SSF
Natural is +0.033 from the Naive Logit, and +0.044 from the SSF
CropsDeveloped Open -0.25 from the Naive Logit, and -0.69 from the SSF

The difference in the interpretation from the different parameters
highlights the ‘effect’ of movement, so to speak, on selection for
covariates. Again, here, note we are not considering any differences
between male or female, or night or day. For example, one could conclude
that the biggest difference is in the dplyr::selection of Developed
Open, which we woudl underestimate the avoidance of if we failed to
consider the movement processes.

## Fitting a ‘Naive’ cLogit Model

Next, we will fit a ‘naive’ clogit model, that is, a model that does not
account for any differences between individuals and treats all
step_id\_’s as independent. Basically ignoring any random effects
structure of individual fishers in this case. \_from
<https://rdrr.io/cran/survival/man/clogit.html_> It turns out that the
loglikelihood for a conditional logistic regression model = loglik from
a Cox model with a particular data structure. Proving this is a nice
homework exercise for a PhD statistics class; not too hard, but the fact
that it is true is surprising.

``` r
clogit1 <- clogit(case_ ~ I(landuse_end) + strata(unique_step), data = fisher6)
#clogit1 <- clogit(case_ ~ I(landuseName=="Developed Open") + I(landuseName=="Developed Other") +I(landuseName=="Natural")+I(landuseName=="Crops") + strata(stratum), data = fisher6)
summary(clogit1)
```

    ## Call:
    ## coxph(formula = Surv(rep(1, 58377L), case_) ~ I(landuse_end) + 
    ##     strata(unique_step), data = fisher6, method = "exact")
    ## 
    ##   n= 58377, number of events= 5307 
    ## 
    ##                     coef exp(coef) se(coef)       z Pr(>|z|)    
    ## I(landuse_end)2 -0.74529   0.47460  0.07760  -9.604  < 2e-16 ***
    ## I(landuse_end)3 -1.75484   0.17293  0.10364 -16.933  < 2e-16 ***
    ## I(landuse_end)5  0.14289   1.15361  0.04572   3.125  0.00178 ** 
    ## I(landuse_end)8 -1.26085   0.28341  0.13691  -9.210  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##                 exp(coef) exp(-coef) lower .95 upper .95
    ## I(landuse_end)2    0.4746     2.1071    0.4076    0.5526
    ## I(landuse_end)3    0.1729     5.7825    0.1411    0.2119
    ## I(landuse_end)5    1.1536     0.8668    1.0547    1.2618
    ## I(landuse_end)8    0.2834     3.5284    0.2167    0.3706
    ## 
    ## Concordance= 0.567  (se = 0.003 )
    ## Likelihood ratio test= 720  on 4 df,   p=<2e-16
    ## Wald test            = 516.1  on 4 df,   p=<2e-16
    ## Score (logrank) test = 584.9  on 4 df,   p=<2e-16

``` r
coef(clogit1)
```

    ## I(landuse_end)2 I(landuse_end)3 I(landuse_end)5 I(landuse_end)8 
    ##      -0.7452918      -1.7548433       0.1428931      -1.2608547

``` r
# tidy up coefficients
clogit_1 <- broom::tidy(clogit1) %>% 
  filter(!term=="(Intercept)") 
## make a figure
figclogit1 <- clogit_1 %>%
  ggplot(., aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = estimate+ 1.96*std.error, ymax = estimate-1.96*std.error)) +
  labs(x = "Habitat", y = "Relative selection Strength") + 
    theme_light() +
  scale_x_discrete(labels = c("Dev(open)", "Dev(other)", "Natural", "Crops")) + geom_hline(yintercept = 0, lty = 2) + ylim(c(-3.75,1.5))
plot_grid(p1, figclogit1)
```

![](README_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

The coefficients are a bit different, but still, quite close even with
the cluster for each individual Stratum.

## Mixed-effect cLogit Models with `coxme`

We can currently fit mixed-effects clogit models using `coxme`. We will
use the `coxme` package to fit a mixed-effects conditional logistic
regression model, accounting for the random effect of individual Fisher
ID in this case. First, we ensure we have the `coxme` package loaded.
See here for more information:
<https://cran.r-project.org/web/packages/coxme/vignettes/coxme.pdf>

``` r
require(coxme)
```

There is no ‘convenient’ wrapper around the `coxme(Surv())` function
like in clogit above. Thus, first, we need to make a ‘fake’ time
variable to trick the Cox-proportional hazards model that time is
irrelevant in your conditional logistic model. We do this by adding a
new variable, time\_ to the fisher6 data frame above. This is actually
what the clogit wrapper around survival is doing, we just don’t know it.

``` r
fisher6$time_ <- ifelse(fisher6$case_ == 0, 2, 1)   #2 for control, 1 for case
table(fisher6$time_, fisher6$case_)
```

    ##    
    ##     FALSE  TRUE
    ##   1     0  5307
    ##   2 53070     0

``` r
head(fisher6)
```

    ## # A tibble: 6 × 18
    ##      id burst_      x1_      x2_    y1_    y2_    sl_    ta_ t1_                
    ##   <dbl>  <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dttm>             
    ## 1  1072      4 1784400. 1784323. 2.41e6 2.41e6  83.3   2.79  2011-02-11 19:50:58
    ## 2  1072      4 1784400. 1784435. 2.41e6 2.41e6  34.8  -0.619 2011-02-11 19:50:58
    ## 3  1072      4 1784400. 1784402. 2.41e6 2.41e6   1.37 -0.942 2011-02-11 19:50:58
    ## 4  1072      4 1784400. 1784438. 2.41e6 2.41e6  75.3  -1.78  2011-02-11 19:50:58
    ## 5  1072      4 1784400. 1784296. 2.41e6 2.41e6 201.   -2.85  2011-02-11 19:50:58
    ## 6  1072      4 1784400. 1784410. 2.41e6 2.41e6 148.   -2.24  2011-02-11 19:50:58
    ## # … with 9 more variables: t2_ <dttm>, dt_ <drtn>, case_ <lgl>, step_id_ <dbl>,
    ## #   landuse_start <dbl>, landuse_end <fct>, unique_step <chr>,
    ## #   landuseName <chr>, time_ <dbl>

``` r
clogitM1<- coxme(Surv(time_,case_) ~ I(landuse_end) + strata(unique_step) + (1|id), data=fisher6)
AIC(clogitM1)
```

    ## [1] 24940.24

``` r
summary(clogitM1)
```

    ## Cox mixed-effects model fit by maximum likelihood
    ##   Data: fisher6
    ##   events, n = 5307, 58377
    ##   Iterations= 5 32 
    ##                     NULL Integrated    Fitted
    ## Log-likelihood -12826.05  -12466.12 -12466.11
    ## 
    ##                    Chisq   df p    AIC    BIC
    ## Integrated loglik 719.86 5.00 0 709.86 676.98
    ##  Penalized loglik 719.87 4.01 0 711.86 685.51
    ## 
    ## Model:  Surv(time_, case_) ~ I(landuse_end) + strata(unique_step) + (1 |      id) 
    ## Fixed coefficients
    ##                       coef exp(coef)   se(coef)      z      p
    ## I(landuse_end)2 -0.7450979 0.4746878 0.07757224  -9.61 0.0000
    ## I(landuse_end)3 -1.7542924 0.1730296 0.10361236 -16.93 0.0000
    ## I(landuse_end)5  0.1423532 1.1529838 0.04565313   3.12 0.0018
    ## I(landuse_end)8 -1.2602773 0.2835754 0.13686634  -9.21 0.0000
    ## 
    ## Random effects
    ##  Group Variable  Std Dev      Variance    
    ##  id    Intercept 8.992044e-03 8.085686e-05

``` r
#clogitM2<- coxme(Surv(time_,case_) ~ I(landuse_end) + strata(stratum) + (1|id) +(I(landuseName=="Developed Other")|id) , data=fisher6)
#summary(clogitM2)
```

We note, importantly, that the random effects variance is quite low,
nearly zero. This means there is very little variation between
individual Fishers in selection. This becomes important later. Moreover
this is very similar to the take home message of Muff et al. (2020) that
using just a random intercept in fitting any RSF or SSF model tends to
underestimate the variance components due to the random intercept or
random coefficients. We will see this below when we fit the same model
with `glmmTMB1`

Second, note, that I tried to fit a similar model using the `mclogit`
package. However, the current implementation of random effects is
limited to the PQL technique, which requires large cluster sizes. Thus,
here, we do not have large enough clusters with only 9 random points. We
see an error message accordingly. Therefore, I have suppressed the
evaluation of the following R code but leave it in for your information
if mcclogit seems useful in the future. The basic call of the mclogit
package requires a column that contains the choice counts or choice
indicators (alternative is chosen=1, is not chosen=0). The second column
contains unique numbers for each choice set. So then the model I tried
to run, but failed, was:

    mclogitTest2 <- mclogit(cbind(used_, unique_step) ~I(landuse_end), random=~ 1|id, data=fisher6)

Plotting the coefficients from `coxme` are more difficult because there
is no tidy approach for objects generated by coxme models, but we can
quickly compare the coefficients directly here, where column 1, 2, and 3
are the naive GLM, naive clogit, and mixed-effect cLogit models. And
compare them to the coefficients for the two-step models from Singer et
al. 

``` r
v1<-model1$coefficients[2:5]
v2<-coef(clogit1)
v3<-coef(clogitM1)
v4<-d2$mean # note that this has female and male two-stage parameter estimates averaged. 
coefSum <- as.data.frame(cbind(v1, v2, v3, v4))
names(coefSum) <- c("Naive", "clogit", "coxme", "two-stage iSSF")
head(coefSum) # 
```

    ##                       Naive     clogit      coxme two-stage iSSF
    ## I(landuse_end)2 -0.42902304 -0.7452918 -0.7450979     -0.9452061
    ## I(landuse_end)3 -1.41033755 -1.7548433 -1.7542924     -2.0644900
    ## I(landuse_end)5  0.05581288  0.1428931  0.1423532      0.1125100
    ## I(landuse_end)8 -1.11393667 -1.2608547 -1.2602773     -1.3217314

Recall that landuse_end 2, 3, 5 and 8 correspond to developed (open),
developed (other), natural, and cropland. And that the forested wet
areas are the ‘reference’ category.

We note that there are some important differences here between the
‘naive’ , clogit/coxme, and two-stage iSSF model. The reason why there
are few differences between clogit, the coxme or mcclogit is because of
the rather low variance in the summary of clogitM1 for individual
Fisher’s we saw above.

However, there are - just like Lab 7 - differences between the
two-staged model, d2, and the coxme/clogit models for landcover types 3
and 5, in particular.

## Model Selection

``` r
AIC(model1, clogit1, clogitM1)
```

    ## Warning in AIC.default(model1, clogit1, clogitM1): models are not all fitted to
    ## the same number of observations

    ##                df      AIC
    ## model1   5.000000 35033.64
    ## clogit1  4.000000 24834.81
    ## clogitM1 4.005735 24940.24

This clearly confirms that really, the mixed-effect model structure is
not necessary and we are fine making inferences using the ‘naive’ cLogit
Model. Finally, we should compare manually the sum of the individual
model AIC values from the two-stage modeling to really undestand if this
clogit model is better than the two-step. However, these results may be
idiosyncratically dependent on the small sample size of Fisher’s here,
6, and the limited variation in response to just one categorical
covariate. Often, addition of mixed-effects models improves model fit
substantially.

## Fitting iSSF models with glmmTMB

Next we will review the Step Selection Function section of Muff et
al. (2020) and follow their advice to fit SSF models using glmmTMB with
a Poisson distribution and random coefficients at the individual level.
This explains in part our noted lack of variance in the random
intercepts for individuals above for the `coxme` models.

Recall: Muff, S., J. Signer, and J. Fieberg. 2020. Accounting for
individual-specific variation in habitat-selection studies: Efficient
estimation of mixed-effects models using Bayesian or frequentist
computation. Journal of Animal Ecology
89:80-92.https://doi.org/10.1111/1365-2656.13087

``` r
head(fisher6)
```

    ## # A tibble: 6 × 18
    ##      id burst_      x1_      x2_    y1_    y2_    sl_    ta_ t1_                
    ##   <dbl>  <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dttm>             
    ## 1  1072      4 1784400. 1784323. 2.41e6 2.41e6  83.3   2.79  2011-02-11 19:50:58
    ## 2  1072      4 1784400. 1784435. 2.41e6 2.41e6  34.8  -0.619 2011-02-11 19:50:58
    ## 3  1072      4 1784400. 1784402. 2.41e6 2.41e6   1.37 -0.942 2011-02-11 19:50:58
    ## 4  1072      4 1784400. 1784438. 2.41e6 2.41e6  75.3  -1.78  2011-02-11 19:50:58
    ## 5  1072      4 1784400. 1784296. 2.41e6 2.41e6 201.   -2.85  2011-02-11 19:50:58
    ## 6  1072      4 1784400. 1784410. 2.41e6 2.41e6 148.   -2.24  2011-02-11 19:50:58
    ## # … with 9 more variables: t2_ <dttm>, dt_ <drtn>, case_ <lgl>, step_id_ <dbl>,
    ## #   landuse_start <dbl>, landuse_end <fct>, unique_step <chr>,
    ## #   landuseName <chr>, time_ <dbl>

``` r
table(fisher6$landuse_end)
```

    ## 
    ##     1     2     3     5     8 
    ## 16259  5411  4877 29972  1858

``` r
TMBm1 <- glmmTMB(case_~ I(landuse_end) + (1|step_id_) + (0 + I(landuse_end)|id), family = poisson, data = fisher6, doFit = FALSE)

TMBm1$parameters$theta[1] <-log(1e3)

TMBm1$mapArg <-list(theta=factor(c(NA, 1:15)))
glmm.TMB.random <- glmmTMB::fitTMB(TMBm1)
summary(glmm.TMB.random)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## case_ ~ I(landuse_end) + (1 | step_id_) + (0 + I(landuse_end) |      id)
    ## Data: fisher6
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##  61872.2  62051.7 -30916.1  61832.2    58357 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups   Name            Variance  Std.Dev.  Corr                    
    ##  step_id_ (Intercept)     1.000e+06 1.000e+03                         
    ##  id       I(landuse_end)1 1.583e-02 1.258e-01                         
    ##           I(landuse_end)2 1.367e-01 3.697e-01  0.19                   
    ##           I(landuse_end)3 1.614e-01 4.017e-01  0.43 -0.08             
    ##           I(landuse_end)5 7.577e-04 2.753e-02  0.30  0.80  0.51       
    ##           I(landuse_end)8 1.522e-01 3.902e-01  0.35  0.61 -0.61  0.09 
    ## Number of obs: 58377, groups:  step_id_, 1788; id, 6
    ## 
    ## Conditional model:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -3.03704   23.64927  -0.128 0.897817    
    ## I(landuse_end)2 -0.58705    0.17759  -3.306 0.000948 ***
    ## I(landuse_end)3 -1.62262    0.19804  -8.193 2.54e-16 ***
    ## I(landuse_end)5  0.07152    0.06195   1.154 0.248343    
    ## I(landuse_end)8 -1.07266    0.20299  -5.284 1.26e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Note that we get an error about model convergence problems with a
non-positive-definitive Hessian Matrix. See here for more details:
<https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html>
This error message is accompanied sometimes by NA’s for SE’s, and
missing/un-estimated LL, AIC, BC and Deviance.

Regardless, we can now add the conditional parameter estimates for the
landcover categories to our list of different model fits we generated
above:

``` r
coefTMB <- fixef(glmm.TMB.random)
v5 <- coefTMB$cond[2:5]
coefSum2 <- as.data.frame(cbind(coefSum, v5))
names(coefSum2) <- c("Naive", "clogit", "coxme", "two-stage iSSF", "glmmTMB")
coefSum2
```

    ##                       Naive     clogit      coxme two-stage iSSF     glmmTMB
    ## I(landuse_end)2 -0.42902304 -0.7452918 -0.7450979     -0.9452061 -0.58705463
    ## I(landuse_end)3 -1.41033755 -1.7548433 -1.7542924     -2.0644900 -1.62262374
    ## I(landuse_end)5  0.05581288  0.1428931  0.1423532      0.1125100  0.07151832
    ## I(landuse_end)8 -1.11393667 -1.2608547 -1.2602773     -1.3217314 -1.07266124

``` r
# recall again that the coefficients are for males and females for two-stage iSSF, which is why there are 8 rows. 
p1
```

![](README_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

In conclusion, there docent seem to be much difference between the
clogit and `coxme` models, but, there is with the `glmmTMB` model.

# SSF Homework

Due April 9 in class.

Conduct an SSF model for wolves in the Cascade, Red Deer, Wildhorse,
Ranch and Bow Valley wolf packs for some covariates that we have used
this semester such as one continuous variable (elevation, slope) and 1
categorical variable (open, closed landcover types). Test for
interactions with step length with covariates. Focus just on one season,
and just 1 level of random effects at the pack level including
intercepts and coefficients as recommended by Muff et al. (2020).
Finally, compare your ssf model estimated with a naive clogit, coxme,
and glmmTMB. Optional: Try your hand at mapping the iSSF using the
movement kernel, transient UD, etc. That might take a long time though.

``` r
wolfGPS <- read.csv("Data/wolfGPS.csv")
head(wolfGPS)
```

    ##     ID1   ID WOLFNAME WOLFUID    PACK FIXID    DATE YEAR SEASON BIOSEASON
    ## 1 23589 3527  Nakomis      85 Cascade  3527 8/22/04 2004 summer       s04
    ## 2 23590 3528  Nakomis      85 Cascade  3528 8/22/04 2004 summer       s04
    ## 3 23591 3529  Nakomis      85 Cascade  3529 8/23/04 2004 summer       s04
    ## 4 23592 3530  Nakomis      85 Cascade  3530 8/23/04 2004 summer       s04
    ## 5 23593 3531  Nakomis      85 Cascade  3531 8/23/04 2004 summer       s04
    ## 6 23594 3532  Nakomis      85 Cascade  3532 8/23/04 2004 summer       s04
    ##       TIME MONTH DAY HOUR JULIANDAY      LAT      LONG ALTITUDE FIXTIME TEMP
    ## 1 18991230     8  22   20       235 51.58880 -115.6574 1836.617      56 15.5
    ## 2 18991230     8  22   22       235 51.59202 -115.6468 1789.390      24  8.5
    ## 3 18991230     8  23    0       236 51.58212 -115.6669 1787.792      59 20.0
    ## 4 18991230     8  23    2       236 51.58183 -115.6665 1815.464      68 22.5
    ## 5 18991230     8  23    4       236 51.58658 -115.6780 1891.054      25  8.0
    ## 6 18991230     8  23    6       236 51.59680 -115.6958 2264.160      53  6.5
    ##             FIXSTATUS N_SATS PDOP X_COORD1 Y_COORD1
    ## 1 "   ""3D Fix-F1"" "      8  1.7 593011.2  5716159
    ## 2 "   ""3D Fix-F1"" "      6  2.3 593738.2  5716531
    ## 3 "   ""3D Fix-F1"" "      5  3.6 592366.6  5715404
    ## 4 "   ""3D Fix-F1"" "      7  5.4 592392.2  5715373
    ## 5 "   ""3D Fix-F1"" "      5  3.9 591585.8  5715886
    ## 6 "   ""3D Fix-F1"" "      6  3.5 590335.7  5717001

``` r
ggplot(wolfGPS, aes(X_COORD1, Y_COORD1, colour = PACK)) +geom_point()
```

![](README_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->
