multiple_testing_adjustment <- function(all, p.adjust.methods, correct, p.values) {
    pvs = unlist(lapply(all, `[[`, "pv"), use.names = FALSE)
    p.rng = !is.na(pvs)
    p = pvs[p.rng]
    lp = length(p)
    min.p = min(p)
    if (correct) {
      pvs.correct = unlist(lapply(all, `[[`, "pv.correct"), use.names = FALSE)
      pv.correct = pvs.correct[p.rng]
      min.pv.correct = min(pv.correct)
    }
    
    if (!is.null(p.adjust.methods)) {
      if ("H" %in% p.adjust.methods) {
        p.values["H"] = min(1, lp * min.p)
      }
      
      if ("Hcorrected" %in% p.adjust.methods) {
        p.values["Hcorrected"] = min(1, lp * min.pv.correct)
      }
    }
    return( p.values = pmin(p.values, 1) )
  }