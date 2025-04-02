load_onearm_data <- function(
    survdat,
    longdat,
    fmla.tte,
    fmla.long,
    longdat.time,
    id.indicator = "id",
    chgpt_location = 0.5, b0_location = 0, b1_location = -0.5, b2_location = 0.5,
    chgpt_scale = 0.5, b0_scale = 1, b1_scale = 0.5, b2_scale = 0.5,
    chgpt_shape = 8, b0_shape = 8, b1_shape = 8, b2_shape = 8) {
  library(data.table)
  library(dplyr)

  # check if survdat and longdat have same "id"
  if (!(id.indicator %in% colnames(survdat) && id.indicator %in% colnames(longdat))) {
    print("The id indicator is not found in one of the datasets.")
    stop()
  }

  tteoutcome.name <- all.vars(fmla.tte)[1]
  tteeventind.name <- all.vars(fmla.tte)[2]
  ttepredictor.name <- all.vars(fmla.tte)[-c(1:2)]
  longpredictor.name <- all.vars(fmla.long)[-1]
  longoutcome.name <- all.vars(fmla.long)[1]

  if (!(tteeventind.name %in% names(longdat))) {
    longdat <- merge(longdat, survdat %>% select(all_of(id.indicator), all_of(tteeventind.name)), by = id.indicator)
  }
  longdat.obs <- longdat %>%
    filter((!!sym(tteeventind.name)) == TRUE, (!!sym(longdat.time)) > 0)
  longdat.obs[["id.new"]] <- as.numeric(as.factor(longdat.obs[[id.indicator]]))
  longdat.obs <- longdat.obs %>%
    arrange(id.new, !!sym(longdat.time))


  longdat.cen <- longdat %>%
    filter((!!sym(tteeventind.name)) == FALSE, (!!sym(longdat.time)) > 0)
  longdat.cen[["id.new"]] <- as.numeric(as.factor(longdat.cen[[id.indicator]]))
  longdat.cen <- longdat.cen %>%
    arrange(id.new, !!sym(longdat.time))

  iddata.obs <- longdat.obs %>%
    group_by(!!sym(id.indicator), id.new) %>%
    summarise(.groups = "drop") %>%
    as.data.frame()
  iddata.cen <- longdat.cen %>%
    group_by(!!sym(id.indicator), id.new) %>%
    summarise(.groups = "drop") %>%
    as.data.frame()

  survdat.obs <- merge(survdat, iddata.obs, by = id.indicator)
  survdat.cen <- merge(survdat, iddata.cen, by = id.indicator)

  standat.list <-
    list(
      nobs = nrow(survdat.obs),
      Nobs = nrow(longdat.obs),
      plong = ncol(model.matrix(fmla.long, data = longdat.cen)),
      ptte = ncol(model.matrix(fmla.tte, data = survdat.cen)),
      idobs = longdat.obs$id.new,
      tobs = survdat.obs[[tteoutcome.name]],
      yobs = c(longdat.obs[[longoutcome.name]]),
      visitobs = longdat.obs[[longdat.time]],
      Xlong_obs = model.matrix(fmla.long, data = longdat.obs),
      Xtte_obs = model.matrix(fmla.tte, data = survdat.obs),
      ncen = nrow(survdat.cen),
      Ncen = nrow(longdat.cen),
      idcen = longdat.cen$id.new,
      tcen = survdat.cen[[tteoutcome.name]],
      ycen = c(longdat.cen[[longoutcome.name]]),
      visitcen = longdat.cen[[longdat.time]],
      Xlong_cen = model.matrix(fmla.long, data = longdat.cen),
      Xtte_cen = model.matrix(fmla.tte, data = survdat.cen),
      randeff_location = c(chgpt_location, b0_location, b1_location, b2_location),
      randeff_scale = c(chgpt_scale, b0_scale, b1_scale, b2_scale),
      randeff_shape = c(chgpt_shape, b0_shape, b1_shape, b2_shape)
    )

  return(standat.list)
}
