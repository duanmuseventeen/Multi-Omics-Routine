require(dplyr)
require(tidyr)

# base--------------------------------------------------------------------------
row(mtcars)
col(mtcars)
diag(as.matrix(mtcars))
#### dplyr----------------------------------------------------------------------
# cheatlist:
arrange(): Arrange rows by variables
distinct(): Select distinct/unique rows
filter(): Return rows with matching conditions
mutate(): Create or transform variables
pull(): Pull out a single variable
relocate(): Change column order
rename(): Rename variables by name
rename_with(): Rename variables with a function
select(): Select variables by name
summarise(): Reduce multiple values down to a single value
slice(): Choose rows by position
# dplyr	base
inner_join(df1, df2)	merge(df1, df2)
left_join(df1, df2)	merge(df1, df2, all.x = TRUE)
right_join(df1, df2)	merge(df1, df2, all.y = TRUE)
full_join(df1, df2)	merge(df1, df2, all = TRUE)
semi_join(df1, df2)	df1[df1$x %in% df2$x, , drop = FALSE]
anti_join(df1, df2)	df1[!df1$x %in% df2$x, , drop = FALSE]
# dplyr::slice------------------------------------------------------------------
mtcars %>% 
  slice(c(1, 2, 5))

slice_head
slice_tail
mtcars %>% slice_max(mpg, n = 5)
mtcars %>% slice_min(mpg, n = 5)
mtcars %>% slice_min(cyl, n = 1, with_ties = FALSE)
mtcars %>% slice_sample
# you can optionally weight by a variable - this code weights by the
# physical weight of the cars, so heavy cars are more likely to get
# selected
mtcars %>% slice_sample(weight_by = wt, n = 5)
# dplyr::column-wise-------------------------------------------------------------------
mtcars %>% 
  group_by(cyl, vs) %>% 
  summarise(mpg = mean(mpg), wt = mean(wt), disp = mean(disp), hp = mean(hp))
  
mtcars %>% 
  group_by(cyl, vs) %>% 
  summarise(across(c('mpg','wt','disp','hp'), mean))
# summarise(across(c(mpg,wt,disp,hp), mean)) # same

mtcars %>% 
  summarise(across(where(is.numeric), n_distinct))


min_max <- list(
  min = ~min(.x, na.rm = TRUE), 
  max = ~max(.x, na.rm = TRUE)
)
starwars %>% summarise(across(where(is.numeric), min_max))
starwars %>% summarise(across(c(height, mass, birth_year), min_max))

starwars %>% summarise(across(where(is.numeric), min_max, .names = "{.fn}.{.col}"))
starwars %>% summarise(across(c(height, mass, birth_year), min_max, .names = "{.fn}.{.col}"))

starwars %>% summarise(
  across(c(height, mass, birth_year), ~min(.x, na.rm = TRUE), .names = "min_{.col}"),
  across(c(height, mass, birth_year), ~max(.x, na.rm = TRUE), .names = "max_{.col}")
)

starwars %>% summarise(
    across(where(is.numeric), ~min(.x, na.rm = TRUE), .names = "min_{.col}"),
    across(where(is.numeric), ~max(.x, na.rm = TRUE), .names = "max_{.col}")  
  )
starwars %>% summarise(
  tibble(
    across(where(is.numeric), ~min(.x, na.rm = TRUE), .names = "min_{.col}"),
    across(where(is.numeric), ~max(.x, na.rm = TRUE), .names = "max_{.col}")  
  )
)

starwars %>% 
  summarise(across(where(is.numeric), min_max, .names = "{.fn}.{.col}")) %>% 
  relocate(starts_with("min"))

df <- tibble(x = 1:3, y = 3:5, z = 5:7)
mult <- list(x = 1, y = 10, z = 100)

df %>% mutate(across(all_of(names(mult)), ~ .x * mult[[cur_column()]]))
#> # A tibble: 3 × 3
#>       x     y     z
#>   <dbl> <dbl> <dbl>
#> 1     1    30   500
#> 2     2    40   600
#> 3     3    50   700

df <- data.frame(x = c(1, 2, 3), y = c(1, 4, 9))
df %>%
  summarise(n = n(), across(where(is.numeric) & !n, sd))

starwars %>% reframe(pick(contains("color")))
starwars %>% distinct(pick(contains("color")))

starwars %>% 
  filter(if_any(everything(), ~ !is.na(.x)))
starwars %>% 
  filter(if_all(everything(), ~ !is.na(.x)))
# dplyr::group-------------------------------------------------------------------------
group_by()
group_keys()
group_indices()
group_rows()
group_vars()
ungroup() 
# dplyr::programming------------------------------------------------------------
# How-tos-----------------------------------------------------------------------
The following examples solve a grab bag of common problems. We show you the minimum amount of code so that you can get the basic idea; most real problems will require more code or combining multiple techniques.

User-supplied data
If you check the documentation, you’ll see that .data never uses data masking or tidy select. That means you don’t need to do anything special in your function:
  
  mutate_y <- function(data) {
    mutate(data, y = a + x)
  }
One or more user-supplied expressions
If you want the user to supply an expression that’s passed onto an argument which uses data masking or tidy select, embrace the argument:
  
  my_summarise <- function(data, group_var) {
    data %>%
      group_by({{ group_var }}) %>%
      summarise(mean = mean(mass))
  }
This generalises in a straightforward way if you want to use one user-supplied expression in multiple places:
  
  my_summarise2 <- function(data, expr) {
    data %>% summarise(
      mean = mean({{ expr }}),
      sum = sum({{ expr }}),
      n = n()
    )
  }
If you want the user to provide multiple expressions, embrace each of them:
  
  my_summarise3 <- function(data, mean_var, sd_var) {
    data %>% 
      summarise(mean = mean({{ mean_var }}), sd = sd({{ sd_var }}))
  }
If you want to use the name of a variable in the output, you can embrace the variable name on the left-hand side of := with {{:
    
    my_summarise4 <- function(data, expr) {
      data %>% summarise(
        "mean_{{expr}}" := mean({{ expr }}),
        "sum_{{expr}}" := sum({{ expr }}),
        "n_{{expr}}" := n()
      )
    }
    my_summarise5 <- function(data, mean_var, sd_var) {
      data %>% 
        summarise(
          "mean_{{mean_var}}" := mean({{ mean_var }}), 
          "sd_{{sd_var}}" := sd({{ sd_var }})
        )
    }
    Any number of user-supplied expressions
    If you want to take an arbitrary number of user supplied expressions, use .... This is most often useful when you want to give the user full control over a single part of the pipeline, like a group_by() or a mutate().
    
    my_summarise <- function(.data, ...) {
      .data %>%
        group_by(...) %>%
        summarise(mass = mean(mass, na.rm = TRUE), height = mean(height, na.rm = TRUE))
    }
    
    starwars %>% my_summarise(homeworld)
    #> # A tibble: 49 × 3
    #>   homeworld    mass height
    #>   <chr>       <dbl>  <dbl>
    #> 1 Alderaan       64   176.
    #> 2 Aleen Minor    15    79 
    #> 3 Bespin         79   175 
    #> 4 Bestine IV    110   180 
    #> # ℹ 45 more rows
    starwars %>% my_summarise(sex, gender)
    #> `summarise()` has grouped output by 'sex'. You can override using the `.groups`
    #> argument.
    #> # A tibble: 6 × 4
    #> # Groups:   sex [5]
    #>   sex            gender      mass height
    #>   <chr>          <chr>      <dbl>  <dbl>
    #> 1 female         feminine    54.7   172.
    #> 2 hermaphroditic masculine 1358     175 
    #> 3 male           masculine   80.2   179.
    #> 4 none           feminine   NaN      96 
    #> # ℹ 2 more rows
    When you use ... in this way, make sure that any other arguments start with . to reduce the chances of argument clashes; see https://design.tidyverse.org/dots-prefix.html for more details.
# Creating multiple columns-----------------------------------------------------
# Sometimes it can be useful for a single expression to return multiple columns. You can do this by returning an unnamed data frame:
  
quantile_df <- function(x, probs = c(0.25, 0.5, 0.75)) {
    tibble(
      val = quantile(x, probs),
      quant = probs
    )
}

# This sort of function is useful inside summarise() and mutate() which allow you to add multiple columns by returning a data frame:
df <- tibble(
    grp = rep(1:3, each = 10),
    x = runif(30),
    y = rnorm(30)
  )

df %>%
  group_by(grp) %>%
  summarise(quantile_df(x, probs = .5))

df %>%
  group_by(grp) %>%
  summarise(across(x:y, ~ quantile_df(.x, probs = .5), .unpack = TRUE))
# Notice that we set .unpack = TRUE inside across(). This tells across() to unpack the data frame returned by quantile_df() into its respective columns, combining the column names of the original columns (x and y) with the column names returned from the function (val and quant).

# If your function returns multiple rows per group, then you’ll need to switch from summarise() to reframe(). summarise() is restricted to returning 1 row summaries per group, but reframe() lifts this restriction:
df %>%
  group_by(grp) %>%
  reframe(across(x:y, quantile_df, .unpack = TRUE))
# Transforming user-supplied variables------------------------------------------
# If you want the user to provide a set of data-variables that are then transformed, use across() and pick():
my_summarise <- function(data, summary_vars) {
    data %>%
      summarise(across({{ summary_vars }}, ~ mean(., na.rm = TRUE)))
  }
starwars %>% 
  group_by(species) %>% 
  my_summarise(c(mass, height))

# You can use this same idea for multiple sets of input data-variables:
my_summarise <- function(data, group_var, summarise_var) {
    data %>%
      group_by(pick({{ group_var }})) %>% 
      summarise(across({{ summarise_var }}, mean))
  }
Use the .names argument to across() to control the names of the output.

my_summarise <- function(data, group_var, summarise_var) {
  data %>%
    group_by(pick({{ group_var }})) %>% 
    summarise(across({{ summarise_var }}, mean, .names = "mean_{.col}"))
}

# Loop over multiple variables--------------------------------------------------
# If you have a character vector of variable names, and want to operate on them with a for loop, index into the special .data pronoun:
  
for (var in names(mtcars)) {
    mtcars %>% count(.data[[var]]) %>% print()
}

# This same technique works with for loop alternatives like the base R apply() family and the purrr map() family:
  
mtcars %>% 
  names() %>% 
  purrr::map(~ count(mtcars, .data[[.x]]))
# (Note that the x in .data[[x]] is always treated as an env-variable; it will never come from the data.)
# dplyr::row-wise---------------------------------------------------------------
df <- tibble(x = 1:2, y = 3:4, z = 5:6)
df %>% rowwise()

# Like group_by(), rowwise() doesn’t really do anything itself; it just changes how the other verbs work. For example, compare the results of mutate() in the following code:
  
df %>% mutate(m = mean(c(x, y, z)))
df %>% rowwise() %>% mutate(m = mean(c(x, y, z)))

df <- tibble(name = c("Mara", "Hadley"), x = 1:2, y = 3:4, z = 5:6)

df %>% 
  rowwise() %>% 
  summarise(m = mean(c(x, y, z)))
# rowwise() is just a special form of grouping, so if you want to remove it from a data frame, just call ungroup().


df %>% 
  rowwise(name) %>% 
  summarise(m = mean(c(x, y, z)))

df <- tibble(id = 1:6, w = 10:15, x = 20:25, y = 30:35, z = 40:45)
rf <- df %>% rowwise(id)
rf %>% mutate(total = sum(c_across(w:z)))

rf %>% 
  mutate(total = sum(c_across(w:z))) %>% 
  ungroup() %>% 
  mutate(across(w:z, ~ . / total))

df %>% mutate(total = rowSums(pick(where(is.numeric), -id)))
# List-columns------------------------------------------------------------------
# rowwise() operations are a natural pairing when you have list-columns. They allow you to avoid explicit loops and/or functions from the apply() or purrr::map() families.

# Motivation--------------------------------------------------------------------
# Imagine you have this data frame, and you want to count the lengths of each element:
df <- tibble(
    x = list(1, 2:3, 4:6)
  )
You might try calling length():
  
# df %>% mutate(l = length(x))
df %>% mutate(l = lengths(x))
# Or if you’re an experienced R programmer, you might know how to apply a function to each element of a list using sapply(), vapply(), or one of the purrr map() functions:
  
df %>% mutate(l = sapply(x, length))
df %>% mutate(l = purrr::map_int(x, length))
# But wouldn’t it be nice if you could just write length(x) and dplyr would figure out that you wanted to compute the length of the element inside of x? Since you’re here, you might already be guessing at the answer: this is just another application of the row-wise pattern.

df %>% 
  rowwise() %>% 
  mutate(l = length(x))
# Subsetting--------------------------------------------------------------------
# Before we continue on, I wanted to briefly mention the magic that makes this work. This isn’t something you’ll generally need to think about (it’ll just work), but it’s useful to know about when something goes wrong.
# There’s an important difference between a grouped data frame where each group happens to have one row, and a row-wise data frame where every group always has one row. Take these two data frames:
  
df <- tibble(g = 1:2, y = list(1:3, "a"))
gf <- df %>% group_by(g)
rf <- df %>% rowwise(g)
# If we compute some properties of y, you’ll notice the results look different:
  
gf %>% mutate(type = typeof(y), length = length(y))
#> # A tibble: 2 × 4
#> # Groups:   g [2]
#>       g y         type  length
#>   <int> <list>    <chr>  <int>
#> 1     1 <int [3]> list       1
#> 2     2 <chr [1]> list       1
rf %>% mutate(type = typeof(y), length = length(y))
#> # A tibble: 2 × 4
#> # Rowwise:  g
#>       g y         type      length
#>   <int> <list>    <chr>      <int>
#> 1     1 <int [3]> integer        3
#> 2     2 <chr [1]> character      1
# They key difference is that when mutate() slices up the columns to pass to length(y) the grouped mutate uses [ and the row-wise mutate uses [[. The following code gives a flavour of the differences if you used a for loop:
# grouped
out1 <- integer(2)
for (i in 1:2) {
  out1[[i]] <- length(df$y[i])
}
out1
#> [1] 1 1

# rowwise
out2 <- integer(2)
for (i in 1:2) {
  out2[[i]] <- length(df$y[[i]])
}
out2
#> [1] 3 1

# Modelling---------------------------------------------------------------------
# rowwise() data frames allow you to solve a variety of modelling problems in what I think is a particularly elegant way. We’ll start by creating a nested data frame:
  
by_cyl <- mtcars %>% nest_by(cyl)
by_cyl
#> # A tibble: 3 × 2
#> # Rowwise:  cyl
#>     cyl data              
#>   <dbl> <list>            
#> 1     4 <tibble [11 × 12]>
#> 2     6 <tibble [7 × 12]> 
#> 3     8 <tibble [14 × 12]>
# This is a little different to the usual group_by() output: we have visibly changed the structure of the data. Now we have three rows (one for each group), and we have a list-col, data, that stores the data for that group. Also note that the output is rowwise(); this is important because it’s going to make working with that list of data frames much easier.

# Once we have one data frame per row, it’s straightforward to make one model per row:
  
mods <- by_cyl %>% mutate(mod = list(lm(mpg ~ wt, data = data)))
mods
#> # A tibble: 3 × 3
#> # Rowwise:  cyl
#>     cyl data               mod   
#>   <dbl> <list>             <list>
#> 1     4 <tibble [11 × 12]> <lm>  
#> 2     6 <tibble [7 × 12]>  <lm>  
#> 3     8 <tibble [14 × 12]> <lm>
# And supplement that with one set of predictions per row:
  
mods <- mods %>% mutate(pred = list(predict(mod, data)))
mods
#> # A tibble: 3 × 4
#> # Rowwise:  cyl
#>     cyl data               mod    pred      
#>   <dbl> <list>             <list> <list>    
#> 1     4 <tibble [11 × 12]> <lm>   <dbl [11]>
#> 2     6 <tibble [7 × 12]>  <lm>   <dbl [7]> 
#> 3     8 <tibble [14 × 12]> <lm>   <dbl [14]>
 
# You could then summarise the model in a variety of ways:
mods %>% summarise(rmse = sqrt(mean((pred - data$mpg) ^ 2)))
#> `summarise()` has grouped output by 'cyl'. You can override using the `.groups`
#> argument.
#> # A tibble: 3 × 2
#> # Groups:   cyl [3]
#>     cyl  rmse
#>   <dbl> <dbl>
#> 1     4 3.01 
#> 2     6 0.985
#> 3     8 1.87
mods %>% summarise(rsq = summary(mod)$r.squared)
#> `summarise()` has grouped output by 'cyl'. You can override using the `.groups`
#> argument.
#> # A tibble: 3 × 2
#> # Groups:   cyl [3]
#>     cyl   rsq
#>   <dbl> <dbl>
#> 1     4 0.509
#> 2     6 0.465
#> 3     8 0.423
mods %>% summarise(broom::glance(mod))
#> `summarise()` has grouped output by 'cyl'. You can override using the `.groups`
#> argument.
#> # A tibble: 3 × 13
#> # Groups:   cyl [3]
#>     cyl r.squared adj.r.squared sigma statistic p.value    df logLik   AIC   BIC
#>   <dbl>     <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>  <dbl> <dbl> <dbl>
#> 1     4     0.509         0.454  3.33      9.32  0.0137     1 -27.7   61.5  62.7
#> 2     6     0.465         0.357  1.17      4.34  0.0918     1  -9.83  25.7  25.5
#> 3     8     0.423         0.375  2.02      8.80  0.0118     1 -28.7   63.3  65.2
#> # ℹ 3 more variables: deviance <dbl>, df.residual <int>, nobs <int>

# Or easily access the parameters of each model:
mods %>% reframe(broom::tidy(mod))
#> # A tibble: 6 × 6
#>     cyl term        estimate std.error statistic    p.value
#>   <dbl> <chr>          <dbl>     <dbl>     <dbl>      <dbl>
#> 1     4 (Intercept)    39.6       4.35      9.10 0.00000777
#> 2     4 wt             -5.65      1.85     -3.05 0.0137    
#> 3     6 (Intercept)    28.4       4.18      6.79 0.00105   
#> 4     6 wt             -2.78      1.33     -2.08 0.0918    
#> # ℹ 2 more rows
                                                                                                                                   
# Repeated function calls-------------------------------------------------------
rowwise() doesn’t just work with functions that return a length-1 vector (aka summary functions); it can work with any function if the result is a list. This means that rowwise() and mutate() provide an elegant way to call a function many times with varying arguments, storing the outputs alongside the inputs.

Simulations
I think this is a particularly elegant way to perform simulations, because it lets you store simulated values along with the parameters that generated them. For example, imagine you have the following data frame that describes the properties of 3 samples from the uniform distribution:
  
  df <- tribble(
    ~ n, ~ min, ~ max,
    1,     0,     1,
    2,    10,   100,
    3,   100,  1000,
  )
You can supply these parameters to runif() by using rowwise() and mutate():
  
  df %>% 
  rowwise() %>% 
  mutate(data = list(runif(n, min, max)))
#> # A tibble: 3 × 4
#> # Rowwise: 
#>       n   min   max data     
#>   <dbl> <dbl> <dbl> <list>   
#> 1     1     0     1 <dbl [1]>
#> 2     2    10   100 <dbl [2]>
#> 3     3   100  1000 <dbl [3]>
Note the use of list() here - runif() returns multiple values and a mutate() expression has to return something of length 1. list() means that we’ll get a list column where each row is a list containing multiple values. If you forget to use list(), dplyr will give you a hint:
  
  df %>% 
  rowwise() %>% 
  mutate(data = runif(n, min, max))
#> Error in `mutate()`:
#> ℹ In argument: `data = runif(n, min, max)`.
#> ℹ In row 2.
#> Caused by error:
#> ! `data` must be size 1, not 2.
#> ℹ Did you mean: `data = list(runif(n, min, max))` ?
# Multiple combinations---------------------------------------------------------
# What if you want to call a function for every combination of inputs? You can use expand.grid() (or tidyr::expand_grid()) to generate the data frame and then repeat the same pattern as above:
  
df <- expand.grid(mean = c(-1, 0, 1), sd = c(1, 10, 100))

df %>% 
  rowwise() %>% 
  mutate(data = list(rnorm(10, mean, sd)))
#> # A tibble: 9 × 3
#> # Rowwise: 
#>    mean    sd data      
#>   <dbl> <dbl> <list>    
#> 1    -1     1 <dbl [10]>
#> 2     0     1 <dbl [10]>
#> 3     1     1 <dbl [10]>
#> 4    -1    10 <dbl [10]>
#> # ℹ 5 more rows
Varying functions
In more complicated problems, you might also want to vary the function being called. This tends to be a bit more of an awkward fit with this approach because the columns in the input tibble will be less regular. But it’s still possible, and it’s a natural place to use do.call():
  
  df <- tribble(
    ~rng,     ~params,
    "runif",  list(n = 10), 
    "rnorm",  list(n = 20),
    "rpois",  list(n = 10, lambda = 5),
  ) %>%
  rowwise()

df %>% 
  mutate(data = list(do.call(rng, params)))
#> # A tibble: 3 × 3
#> # Rowwise: 
#>   rng   params           data      
#>   <chr> <list>           <list>    
#> 1 runif <named list [1]> <dbl [10]>
#> 2 rnorm <named list [1]> <dbl [20]>
#> 3 rpois <named list [2]> <int [10]>
Previously
rowwise()
rowwise() was also questioning for quite some time, partly because I didn’t appreciate how many people needed the native ability to compute summaries across multiple variables for each row. As an alternative, we recommended performing row-wise operations with the purrr map() functions. However, this was challenging because you needed to pick a map function based on the number of arguments that were varying and the type of result, which required quite some knowledge of purrr functions.

# I was also resistant to rowwise() because I felt like automatically switching between [ to [[ was too magical in the same way that automatically list()-ing results made do() too magical. I’ve now persuaded myself that the row-wise magic is good magic partly because most people find the distinction between [ and [[ mystifying and rowwise() means that you don’t need to think about it.
                                                                                                                                                                                                                                                                                                                            
Since rowwise() clearly is useful it is not longer questioning, and we expect it to be around for the long term.
do()
We’ve questioned the need for do() for quite some time, because it never felt very similar to the other dplyr verbs. It had two main modes of operation:
  
  Without argument names: you could call functions that input and output data frames using . to refer to the “current” group. For example, the following code gets the first row of each group:
  
  mtcars %>% 
  group_by(cyl) %>% 
  do(head(., 1))
#> # A tibble: 3 × 13
#> # Groups:   cyl [3]
#>     mpg   cyl  disp    hp  drat    wt  qsec    vs    am  gear  carb  cyl2  cyl4
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1  22.8     4   108    93  3.85  2.32  18.6     1     1     4     1     8    16
#> 2  21       6   160   110  3.9   2.62  16.5     0     1     4     4    12    24
#> 3  18.7     8   360   175  3.15  3.44  17.0     0     0     3     2    16    32
This has been superseded by pick() plus reframe(), a variant of summarise() that can create multiple rows and columns per group.

mtcars %>% 
  group_by(cyl) %>% 
  reframe(head(pick(everything()), 1))
#> # A tibble: 3 × 13
#>     cyl   mpg  disp    hp  drat    wt  qsec    vs    am  gear  carb  cyl2  cyl4
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1     4  22.8   108    93  3.85  2.32  18.6     1     1     4     1     8    16
#> 2     6  21     160   110  3.9   2.62  16.5     0     1     4     4    12    24
#> 3     8  18.7   360   175  3.15  3.44  17.0     0     0     3     2    16    32
With arguments: it worked like mutate() but automatically wrapped every element in a list:
  
  mtcars %>% 
  group_by(cyl) %>% 
  do(nrows = nrow(.))
#> # A tibble: 3 × 2
#> # Rowwise: 
#>     cyl nrows    
#>   <dbl> <list>   
#> 1     4 <int [1]>
#> 2     6 <int [1]>
#> 3     8 <int [1]>
I now believe that behaviour is both too magical and not very useful, and it can be replaced by summarise() and pick().

mtcars %>% 
  group_by(cyl) %>% 
  summarise(nrows = nrow(pick(everything())))
#> # A tibble: 3 × 2
#>     cyl nrows
#>   <dbl> <int>
#> 1     4    11
#> 2     6     7
#> 3     8    14
If needed (unlike here), you can wrap the results in a list yourself.

The addition of pick()/across() and the increased scope of summarise()/reframe() means that do() is no longer needed, so it is now superseded.
# dplyr::Two-table verbs--------------------------------------------------------
flights2 %>% left_join(airports, c("dest" = "faa"))
# dplyr::window functions-------------------------------------------------------
# tidyr::separate---------------------------------------------------------------
df <- mtcars %>% 
  mutate(myvar = "A-B-C",
         mynum = "1-2")

df %>% separate(myvar, c("A", "B"), remove = FALSE) # C is not split
df %>% separate(myvar, c("A", "B", "C"), remove = FALSE)
df %>% separate(myvar, c("A", "B", "C", "D"), remove = FALSE) # NA is introduced
df %>% separate(mynum, c("A", "B"), remove = FALSE, convert = TRUE)

df %>% separate_wider_delim(myvar, 
                            delim = "-", 
                            names = c("A", "B", "C"), 
                            cols_remove = FALSE)

df %>% separate_wider_delim(
  myvar, delim = "-", names = c("A", "B"), cols_remove = FALSE, too_many = "drop")
df %>% separate_wider_delim(
  myvar, delim = "-", names = c("A", "B"), cols_remove = FALSE, too_many = "merge")

df %>% separate_wider_delim(
  myvar, delim = "-", names = c("A", "B", "C", "D"), cols_remove = FALSE, too_few = "align_start")
df %>% separate_wider_delim(
  myvar, delim = "-", names = c("A", "B", "C", "D"), cols_remove = FALSE, too_few = "align_end")
# purrr::map--------------------------------------------------------------------
require(purrr)

mtcars |>
  split(mtcars$cyl) |>
  map(\(df) lm(mpg ~ wt, data = df)) |>
  map_dfr(\(mod) as.data.frame(t(as.matrix(coef(mod)))))








