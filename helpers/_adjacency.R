# _adjacency.R
# Graph utilities, wormhole definitions, and build_augmented_adjacency()
# Sourced by 09_fit_model.Rmd

library(RANN)
library(igraph)

# ═══════════════════════════════════════════════════════════════════════════════
# Wormhole Definitions
# ═══════════════════════════════════════════════════════════════════════════════

NYC_WORMHOLES <- list(
  # L0 critical: airports, ferries
  list(name="LGA-Midtown",           from=c(-73.8740,40.7769), to=c(-73.9857,40.7484), priority="critical"),
  list(name="LGA-EastElmhurst",      from=c(-73.8740,40.7769), to=c(-73.8700,40.7610), priority="critical"),
  list(name="JFK-DowntownBK",         from=c(-73.7781,40.6413), to=c(-73.9857,40.6892), priority="critical"),
  list(name="JFK-Jamaica",            from=c(-73.7781,40.6413), to=c(-73.7930,40.7028), priority="critical"),
  list(name="JFK-HowardBeach",        from=c(-73.7781,40.6413), to=c(-73.8300,40.6590), priority="critical"),
  list(name="StGeorge-Whitehall",     from=c(-74.0766,40.6437), to=c(-74.0130,40.7013), priority="critical"),
  list(name="Greenpt-E34thFerry",     from=c(-73.9615,40.7320), to=c(-73.9718,40.7440), priority="critical"),
  # L3 manual: rail, highways
  list(name="BKCruiseTerminal-DTBK",  from=c(-74.0168,40.6808), to=c(-73.9857,40.6892), priority="manual"),
  list(name="BKCruiseTerminal-BayR",  from=c(-74.0168,40.6808), to=c(-74.0302,40.6346), priority="manual"),
  list(name="SunnysideYd-PennSta",    from=c(-73.9215,40.7425), to=c(-73.9937,40.7506), priority="manual"),
  list(name="AtlanticYd-Atlantic",    from=c(-73.9762,40.6846), to=c(-73.9776,40.6862), priority="manual"),
  list(name="GrandCentral-Harlem125", from=c(-73.9772,40.7527), to=c(-73.9390,40.8043), priority="manual"),
  list(name="Fulton-Barclays",        from=c(-74.0071,40.7092), to=c(-73.9762,40.6846), priority="manual"),
  list(name="CrossBronx-GWBridge",    from=c(-73.8800,40.8380), to=c(-73.9520,40.8510), priority="manual"),
  list(name="BQE-HughCarey",          from=c(-74.0020,40.6870), to=c(-74.0145,40.6895), priority="manual")
)

# Utah wormholes: interstate highway junctions connecting non-adjacent counties
UT_WORMHOLES <- list(
  # L0 critical: SLC Airport -> major corridors
  list(name="SLCairport-DowntownSLC", from=c(-111.9791,40.7899), to=c(-111.8910,40.7608), priority="critical"),
  list(name="SLCairport-WestValley",  from=c(-111.9791,40.7899), to=c(-111.9390,40.6916), priority="critical"),
  # I-15 long-distance corridor
  list(name="I15-Weber-Davis",        from=c(-111.9730,41.2230), to=c(-111.9020,40.9880), priority="critical"),
  list(name="I15-Utah-Juab",          from=c(-111.6580,40.2330), to=c(-111.7480,39.7310), priority="critical"),
  list(name="I15-Iron-Washington",    from=c(-113.0610,37.6770), to=c(-113.5080,37.1080), priority="critical"),
  # I-80 corridor
  list(name="I80-SL-Tooele",          from=c(-111.8910,40.7608), to=c(-112.3260,40.6220), priority="critical"),
  list(name="I80-Summit-Wasatch",     from=c(-111.4980,40.7420), to=c(-111.5010,40.6540), priority="critical"),
  # L3 manual: secondary corridors
  list(name="I70-Sevier-Emery",       from=c(-111.8240,38.7830), to=c(-111.2540,38.9960), priority="manual"),
  list(name="US89-Sanpete-Sevier",    from=c(-111.5790,39.3640), to=c(-111.8240,38.7830), priority="manual"),
  list(name="US6-Carbon-Duchesne",    from=c(-110.9530,39.6070), to=c(-110.4030,40.1630), priority="manual"),
  list(name="I84-Weber-BoxElder",     from=c(-111.9730,41.2230), to=c(-112.0150,41.7390), priority="manual")
)


# ═══════════════════════════════════════════════════════════════════════════════
# Graph Utilities
# ═══════════════════════════════════════════════════════════════════════════════

# Iterative component finder (avoids spdep recursion stack overflow on >15k nodes)
safe_n_comp_nb <- function(nb_obj) {
  n <- length(nb_obj)
  el <- vector("list", n)
  for (i in seq_len(n)) {
    nbrs <- nb_obj[[i]]
    nbrs <- nbrs[nbrs > 0L & nbrs <= n]
    if (length(nbrs) > 0) el[[i]] <- data.frame(from = i, to = nbrs)
  }
  edges <- do.call(rbind, el)
  g <- make_empty_graph(n = n, directed = FALSE)
  if (!is.null(edges) && nrow(edges) > 0)
    g <- add_edges(g, as.vector(t(edges[, 1:2])))
  mem <- components(g)
  list(nc = mem$no, comp.id = as.integer(mem$membership))
}

# Symmetric edge insert with 0L island-marker handling
nb_set_edge <- function(nb, i, j) {
  ni <- nb[[i]]; ni <- ni[ni > 0L]
  nj <- nb[[j]]; nj <- nj[nj > 0L]
  nb[[i]] <- sort(unique(c(ni, as.integer(j))))
  nb[[j]] <- sort(unique(c(nj, as.integer(i))))
  nb
}

sanitize_nb <- function(nb_obj) {
  n <- length(nb_obj)
  for (i in seq_len(n)) {
    nbrs <- nb_obj[[i]]
    nbrs <- as.integer(nbrs[!is.na(nbrs)])
    nbrs <- sort(unique(nbrs[nbrs != i & nbrs > 0L & nbrs <= n]))
    nb_obj[[i]] <- if (length(nbrs) == 0) 0L else nbrs
  }
  attr(nb_obj, "class") <- "nb"
  attr(nb_obj, "region.id") <- as.character(1:n)
  nb_obj
}

safe_nb2INLA <- function(graph_file, nb_obj) {
  nb2INLA(graph_file, sanitize_nb(nb_obj))
}


# ═══════════════════════════════════════════════════════════════════════════════
# Wormhole Application Helper
# ═══════════════════════════════════════════════════════════════════════════════

apply_wormholes <- function(nb, wh_list, coords_proj, layer_name) {
  n_wh <- 0
  if (length(wh_list) == 0) return(list(nb = nb, n = 0))
  comp <- safe_n_comp_nb(nb)
  for (wh in wh_list) {
    if (n_wh >= MAX_WORMHOLES) break
    i1 <- which.min(sqrt((coords_proj[,1] - wh$from[1])^2 + (coords_proj[,2] - wh$from[2])^2))
    i2 <- which.min(sqrt((coords_proj[,1] - wh$to[1])^2   + (coords_proj[,2] - wh$to[2])^2))
    if (i1 != i2 && comp$comp.id[i1] != comp$comp.id[i2] && !(i2 %in% nb[[i1]])) {
      nb <- nb_set_edge(nb, i1, i2)
      n_wh <- n_wh + 1
      comp <- safe_n_comp_nb(nb)
    }
  }
  list(nb = nb, n = n_wh)
}


# ═══════════════════════════════════════════════════════════════════════════════
# Augmented Adjacency Builder (10-stage pipeline)
#
# L0 queen -> L1 critical WH -> L2 bridges -> L3 snap ->
# L4 rail -> L5 roads -> L6 OD -> L7 manual WH -> absorption sweep -> classify
# ═══════════════════════════════════════════════════════════════════════════════

build_augmented_adjacency <- function(grid_sf, bridge_links, wormholes, gt_dt = NULL) {
  n_grid <- nrow(grid_sf)
  grid_ids <- grid_sf$grid_id
  utm_crs <- st_crs(grid_sf)

  centroids_proj <- st_centroid(grid_sf)
  coords_proj <- st_coordinates(centroids_proj)

  # Convert wormhole coordinates from WGS84 to UTM
  wh_pts <- do.call(rbind, lapply(wormholes, function(w) rbind(w$from, w$to)))
  wh_sf <- st_as_sf(data.frame(x = wh_pts[,1], y = wh_pts[,2]),
                    coords = c("x", "y"), crs = 4326) |> st_transform(utm_crs)
  wh_utm <- st_coordinates(wh_sf)
  for (i in seq_along(wormholes)) {
    wormholes[[i]]$from <- wh_utm[2*i - 1, 1:2]
    wormholes[[i]]$to   <- wh_utm[2*i, 1:2]
  }

  # Convert bridge coordinates from WGS84 to UTM
  if (nrow(bridge_links) > 0) {
    br_sf <- st_as_sf(bridge_links[!is.na(lat) & !is.na(lon)],
                      coords = c("lon", "lat"), crs = 4326) |> st_transform(utm_crs)
    br_utm <- st_coordinates(br_sf)
    bridge_links[!is.na(lat) & !is.na(lon), c("utm_x", "utm_y") := .(br_utm[,1], br_utm[,2])]
  }

  # ── Stage L0: Queen contiguity ──
  suppressWarnings({ nb <- poly2nb(grid_sf, queen = TRUE, snap = SNAP_QUEEN) })
  comp0 <- safe_n_comp_nb(nb)

  # ── Stage L1: Critical wormholes (airports, ferries) ──
  wh_critical <- Filter(function(w) w$priority == "critical", wormholes)
  wh_manual   <- Filter(function(w) w$priority != "critical", wormholes)
  res0 <- apply_wormholes(nb, wh_critical, coords_proj, "L0")
  nb <- res0$nb

  # ── Stage L2: Bridge & tunnel augmentation ──
  n_bridge <- 0
  if (nrow(bridge_links) > 0) {
    for (i in 1:nrow(bridge_links)) {
      bx <- bridge_links$utm_x[i]; by <- bridge_links$utm_y[i]
      if (is.na(bx) || is.na(by)) next
      dists <- sqrt((coords_proj[,1] - bx)^2 + (coords_proj[,2] - by)^2)
      nearest_idx <- order(dists)[1:min(20, n_grid)]
      cur <- safe_n_comp_nb(nb)
      nc <- cur$comp.id[nearest_idx]
      for (j in 1:(length(nearest_idx) - 1)) {
        linked <- FALSE
        for (k in (j + 1):min(length(nearest_idx), j + 10)) {
          if (nc[j] != nc[k]) {
            ia <- nearest_idx[j]; ib <- nearest_idx[k]
            if (!(ib %in% nb[[ia]])) {
              nb <- nb_set_edge(nb, ia, ib); n_bridge <- n_bridge + 1
            }
            linked <- TRUE; break
          }
        }
        if (linked) break
      }
    }
  }
  # ── Stage L3: Progressive distance snapping ──
  n_snap <- 0
  for (ti in seq_along(SNAP_TIERS)) {
    snap_dist <- min(SNAP_TIERS[ti], SNAP_MAX)
    cur <- safe_n_comp_nb(nb)
    if (cur$nc <= 1) break
    cs <- table(cur$comp.id)
    lg_id <- as.integer(names(which.max(cs)))
    main_g <- which(cur$comp.id == lg_id)
    small_ids <- as.integer(names(cs[cs < 50 & names(cs) != as.character(lg_id)]))
    if (length(small_ids) == 0) break
    n_this <- 0
    for (cid in small_ids) {
      frag <- which(cur$comp.id == cid)
      frag_s <- frag[seq(1, length(frag), length.out = min(10, length(frag)))]
      nn <- nn2(coords_proj[main_g, , drop = FALSE], coords_proj[frag_s, , drop = FALSE], k = 1)
      bi <- which.min(nn$nn.dists)
      if (nn$nn.dists[bi] <= snap_dist) {
        ia <- frag_s[bi]; ib <- main_g[nn$nn.idx[bi]]
        if (!(ib %in% nb[[ia]])) { nb <- nb_set_edge(nb, ia, ib); n_snap <- n_snap + 1; n_this <- n_this + 1 }
      }
    }
    if (n_this == 0 && ti > 2) break
  }

  # ── Stage L4: Rail/subway corridors ──
  n_rail <- 0
  subway_loaded <- FALSE
  if (file.exists(SUBWAY_ROUTES_SHP)) {
    tryCatch({
      subway <- st_transform(st_read(SUBWAY_ROUTES_SHP, quiet = TRUE), st_crs(grid_sf))
      hits <- st_intersects(grid_sf, st_buffer(subway, 75), sparse = TRUE)
      cur <- safe_n_comp_nb(nb)
      for (li in 1:nrow(subway)) {
        touching <- which(lengths(lapply(hits, function(x) x[x == li])) > 0)
        if (length(touching) < 2) next
        lc <- cur$comp.id[touching]; uc <- unique(lc)
        if (length(uc) <= 1) next
        for (ci in 2:length(uc)) {
          g1 <- touching[which(lc == uc[1])[1]]
          g2 <- touching[which(lc == uc[ci])[1]]
          if (!(g2 %in% nb[[g1]])) { nb <- nb_set_edge(nb, g1, g2); n_rail <- n_rail + 1 }
        }
        if (li %% 50 == 0) cur <- safe_n_comp_nb(nb)
      }
      subway_loaded <- TRUE
    }, error = function(e) NULL)
  }
  if (!subway_loaded && exists("SUBWAY_STATIONS_CSV") && file.exists(SUBWAY_STATIONS_CSV)) {
    tryCatch({
      mta <- fread(SUBWAY_STATIONS_CSV)
      mta <- mta[!is.na(`GTFS Latitude`) & !is.na(`GTFS Longitude`)]
      mta_sf <- st_transform(st_as_sf(mta, coords = c("GTFS Longitude", "GTFS Latitude"), crs = 4326), st_crs(grid_sf))
      mta_coords <- st_coordinates(mta_sf)
      cur <- safe_n_comp_nb(nb)
      for (ln in unique(mta$Line)) {
        ln_idx <- which(mta$Line == ln)
        if (length(ln_idx) < 2) next
        grids_on_line <- vapply(ln_idx, function(si) {
          which.min(sqrt((coords_proj[,1] - mta_coords[si,1])^2 + (coords_proj[,2] - mta_coords[si,2])^2))
        }, integer(1))
        for (si in 1:(length(grids_on_line) - 1)) {
          ia <- grids_on_line[si]; ib <- grids_on_line[si + 1]
          if (ia != ib && cur$comp.id[ia] != cur$comp.id[ib]) {
            nb <- nb_set_edge(nb, ia, ib); n_rail <- n_rail + 1
          }
        }
      }
    }, error = function(e) NULL)
  }

  # ── Stage L5: Road network augmentation ──
  n_road <- 0
  road_files <- if (file.exists(ROAD_SHP)) ROAD_SHP else get_road_paths(state)
  for (rf in road_files) {
    if (!file.exists(rf)) next
    tryCatch({
      roads_sf <- st_transform(st_read(rf, quiet = TRUE), st_crs(grid_sf))
      mtfcc_col <- intersect(c("MTFCC", "mtfcc"), names(roads_sf))[1]
      if (!is.na(mtfcc_col)) roads_major <- roads_sf[roads_sf[[mtfcc_col]] %in% TIGER_BRIDGE_TYPES, ]
      else roads_major <- roads_sf
      comp_pre <- safe_n_comp_nb(nb)
      if (comp_pre$nc > 1 && nrow(roads_major) > 0) {
        road_hits <- st_intersects(roads_major, grid_sf)
        for (r in seq_along(road_hits)) {
          tidx <- road_hits[[r]]
          if (length(tidx) < 2) next
          tcomps <- comp_pre$comp.id[tidx]; ucomps <- unique(tcomps)
          if (length(ucomps) < 2) next
          for (ci in 1:(length(ucomps) - 1)) {
            g_c1 <- tidx[tcomps == ucomps[ci]]; g_c2 <- tidx[tcomps == ucomps[ci + 1]]
            best_d <- Inf; best_i <- NA; best_j <- NA
            for (gi in g_c1[1:min(5, length(g_c1))]) {
              for (gj in g_c2[1:min(5, length(g_c2))]) {
                d <- sqrt((coords_proj[gi,1] - coords_proj[gj,1])^2 + (coords_proj[gi,2] - coords_proj[gj,2])^2)
                if (d < best_d) { best_d <- d; best_i <- gi; best_j <- gj }
              }
            }
            if (!is.na(best_i) && best_i != best_j && !(best_j %in% nb[[best_i]])) {
              nb <- nb_set_edge(nb, best_i, best_j); n_road <- n_road + 1
              comp_pre$comp.id[comp_pre$comp.id == ucomps[ci + 1]] <- ucomps[ci]
            }
          }
        }
      }
    }, error = function(e) NULL)
  }

  # ── Stage L6: OD correlation (isolates only) ──
  n_od <- 0
  comp2 <- safe_n_comp_nb(nb)
  if (!is.null(gt_dt) && nrow(gt_dt) > OD_MIN_ROWS && comp2$nc > 1) {
    cs <- table(comp2$comp.id)
    lg_id <- as.integer(names(which.max(cs)))
    main_idx <- which(comp2$comp.id == lg_id)
    main_gids <- grid_ids[main_idx]

    # identify structural isolates
    isolate_grids <- integer(0)
    for (fc in unique(comp2$comp.id[comp2$comp.id != lg_id])) {
      fc_idx <- which(comp2$comp.id == fc)
      fc_deg <- mean(sapply(fc_idx, function(i) length(nb[[i]][nb[[i]] > 0L])))
      if (length(fc_idx) <= OD_MAX_FRAG_SIZE || fc_deg < OD_MAX_FRAG_DEGREE)
        isolate_grids <- c(isolate_grids, fc_idx)
    }
    isolate_gids <- grid_ids[isolate_grids]

    rho_col <- if ("log_rho" %in% names(gt_dt)) "log_rho" else "Y"
    act_col <- if ("rho_norm" %in% names(gt_dt)) "rho_norm" else "Y"
    gstats <- gt_dt[, .(total_act = sum(get(act_col), na.rm = TRUE),
                         n_days = sum(get(act_col) > 0, na.rm = TRUE)), by = grid_id]
    elig <- gstats[grid_id %in% isolate_gids & total_act >= OD_MIN_ACTIVITY & n_days >= OD_MIN_DAYS, grid_id]

    if (length(elig) > 0) {
      fp <- dcast(gt_dt[grid_id %in% elig], grid_id ~ time_id, value.var = rho_col, fill = 0, fun.aggregate = mean)
      fg <- fp$grid_id; fm <- as.matrix(fp[, -1, with = FALSE])
      ms <- gstats[grid_id %in% main_gids][order(-total_act)]
      ms_ids <- ms$grid_id[1:min(300, nrow(ms))]
      mp <- dcast(gt_dt[grid_id %in% ms_ids], grid_id ~ time_id, value.var = rho_col, fill = 0, fun.aggregate = mean)
      mg <- mp$grid_id; mm <- as.matrix(mp[, -1, with = FALSE])
      for (fi in 1:nrow(fm)) {
        fgid <- fg[fi]; fidx <- which(grid_ids == fgid)
        if (length(fidx) == 0) next
        fv <- fm[fi, ]; if (sd(fv, na.rm = TRUE) < 1e-10) next
        cors <- apply(mm, 1, function(mv) if (sd(mv, na.rm = TRUE) < 1e-10) 0 else cor(fv, mv, use = "complete.obs"))
        bi <- which.max(cors)
        if (!is.na(cors[bi]) && cors[bi] > OD_CORR_THRESHOLD) {
          midx <- which(grid_ids == mg[bi])
          if (length(midx) > 0 && !(midx %in% nb[[fidx]])) { nb <- nb_set_edge(nb, fidx, midx); n_od <- n_od + 1 }
        }
      }
    }
  }
  # ── Stage L7: Manual wormholes (only if still fragmented) ──
  n_wh_l3 <- 0
  if (length(wh_manual) > 0 && safe_n_comp_nb(nb)$nc > 1) {
    res3 <- apply_wormholes(nb, wh_manual, coords_proj, "L3")
    nb <- res3$nb; n_wh_l3 <- res3$n
  }
  # ── Stage ABS: Absorption sweep ──
  n_abs <- 0
  for (pass in 1:10) {
    sc <- safe_n_comp_nb(nb)
    if (sc$nc <= 1) break
    cs <- table(sc$comp.id)
    lg_id <- as.integer(names(which.max(cs)))
    main_g <- which(sc$comp.id == lg_id)
    n_this <- 0
    for (cid in unique(sc$comp.id[sc$comp.id != lg_id])) {
      frag <- which(sc$comp.id == cid)
      frag_s <- frag[seq(1, length(frag), length.out = min(5, length(frag)))]
      nn <- nn2(coords_proj[main_g, , drop = FALSE], coords_proj[frag_s, , drop = FALSE], k = 1)
      bi <- which.min(nn$nn.dists)
      ia <- frag_s[bi]; ib <- main_g[nn$nn.idx[bi]]
      if (!(ib %in% nb[[ia]])) { nb <- nb_set_edge(nb, ia, ib); n_abs <- n_abs + 1; n_this <- n_this + 1 }
    }
    if (n_this == 0) break
  }

  # ── Island fixes ──
  n_isl <- 0
  zeros <- which(card(nb) == 0)
  if (length(zeros) > 0) {
    others <- setdiff(1:n_grid, zeros)
    for (iso in zeros) {
      nn <- nn2(coords_proj[others, , drop = FALSE], matrix(coords_proj[iso, ], nrow = 1), k = 1)
      nb <- nb_set_edge(nb, iso, others[nn$nn.idx[1]]); n_isl <- n_isl + 1
    }
  }
  comp_f <- safe_n_comp_nb(nb)
  list(nb = nb, comp = comp_f, comp_before = comp0$nc, comp_after = comp_f$nc,
       n_wh_l0 = res0$n, n_bridge = n_bridge, n_snap = n_snap, n_rail = n_rail,
       n_road = n_road, n_od = n_od, n_wormhole = res0$n + n_wh_l3,
       n_absorbed = n_abs, n_island_fix = n_isl)
}
