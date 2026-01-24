# glyparse (0.5.3)

* GitHub: <https://github.com/glycoverse/glyparse>
* Email: <mailto:23110220018@m.fudan.edu.cn>
* GitHub mirror: <https://github.com/cran/glyparse>

Run `revdepcheck::revdep_details(, "glyparse")` for more info

## Newly broken

*   checking tests ...
     ```
       Running ‘testthat.R’
      ERROR
     Running the tests in ‘tests/testthat.R’ failed.
     Last 13 lines of output:
         2.   └─glyparse:::struc_parser_wrapper(x, do_auto_parse)
         3.     └─glyrepr::as_glycan_structure(unique_graphs)
         4.       └─vctrs::vec_cast(x, glycan_structure())
         5.         └─vctrs (local) `<fn>`()
         6.           └─glyrepr:::vec_cast.glyrepr_structure.list(...)
         7.             ├─base::do.call(glycan_structure, x)
         8.             └─glyrepr (local) `<fn>`(`<igraph>`, `<igraph>`, `<igraph>`)
         9.               └─glyrepr:::validate_glycan_structure_vector(reordered_graphs)
        10.                 └─cli::cli_abort(...)
        11.                   └─rlang::abort(...)
       
       [ FAIL 1 | WARN 0 | SKIP 12 | PASS 291 ]
       Error:
       ! Test failures.
       Execution halted
     ```

*   checking running R code from vignettes ...
     ```
       ‘glyparse.Rmd’ using ‘UTF-8’... failed
      ERROR
     Errors in running code in vignettes:
     when running code in ‘glyparse.Rmd’
       ...
     > x <- c("Gal(b1-3)GalNAc(b1-", "(N(F)(N(H(H(N))(H(N(H))))))", 
     +     "WURCS=2.0/3,3,2/[a2122h-1b_1-5][a1122h-1b_1-5][a1122h-1a_1-5]/1-2-3/a4-b1_b3-c1 ..." ... [TRUNCATED] 
     
     > auto_parse(x)
     
       When sourcing ‘glyparse.R’:
     错误: All structures must have the same monosaccharide type.
     ✖ Found 2 concrete and 1 generic structure(s) in the same vector.
     ℹ Use `convert_to_generic()` to convert concrete structures to generic type.
     停止执行
     ```

