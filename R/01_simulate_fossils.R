############################################################
## 01_simulate_fossils.R
##
## WHAT THIS SCRIPT DOES
## ----------------------------------------------------------
## Runs the full simulation study: builds birth-death trees,
## evolves a binary character on them, generates a
## state-dependent fossil record, reconstructs ancestral
## states from the fossil-inclusive tree, and scores how
## accurate that reconstruction is against the known truth.
##
## It repeats this 500 times for each of 16 preservation
## scenarios (5 heterogeneity levels x 3 overall preservation
## levels, plus a no-fossil control) = 8,000 simulations total.
##
## WARNING: this takes a long time to run (likely
## several hours depending on your machine), because it is
## simulating and analysing 8,000 individual phylogenies.
## You do NOT need to re-run this to explore the results -- see
## 02_analyze_results.R, which works directly from the saved
## CSV output.
##
## OUTPUT
## ----------------------------------------------------------
## simulation_results_revised.csv -- one row per successful
## simulation, saved incrementally every 100 rows (in case the
## run is interrupted) and once more at the end.
############################################################

## ---- Packages ---------------------------------------------
# install.packages("phangorn")  # run once if not already installed

library(phytools)
library(phangorn)
library(TreeSim)
library(FossilSim)
library(paleotree)

# `treesurgeon` was loaded in the original script but is not
# actually used anywhere below. Leaving it commented out rather
# than deleting the line outright -- uncomment if you find it's
# needed for something downstream that isn't obvious from this
# script alone.
# library(treesurgeon)


############################################################
## SECTION 1: Helper functions
##
## Internal utility functions used to construct the simulated
## fossil tree. You shouldn't need to call these directly --
## they're used internally by simulate_fossils() and
## build_fossil_tree() below.
############################################################

## Splits a single simmap branch (a named vector of state
## durations) into two pieces at a given distance from the
## start of the branch. Used when a fossil falls partway along
## a branch and we need to know the state at that exact point.
split_map <- function(map, distance) {
  cum <- cumsum(map)

  i <- which(cum >= distance)[1]

  before <- map[1:i]
  after <- map[i:length(map)]

  ## Split the interval containing the split point
  excess <- cum[i] - distance

  before[length(before)] <- before[length(before)] - excess
  after[1] <- excess

  ## Remove zero-length segments created by rounding
  before <- before[before > 1e-12]
  after <- after[after > 1e-12]

  list(
    upper = before,
    lower = after,
    state = names(before)[length(before)]
  )
}

## Identifies the chain of edges in the fossil tree that
## correspond to a single original branch, after fossils have
## been inserted along it as zero-length tips.
get_edge_chain <- function(f_tree, tip_names) {
  ## Starting node
  if (length(tip_names) == 1) {
    current_node <- which(f_tree$tip.label == tip_names)
  } else {
    current_node <- mrca.phylo(f_tree, tip_names)
  }

  edge_chain <- integer()

  repeat {
    ## Incoming edge
    edge <- which(f_tree$edge[, 2] == current_node)

    ## Stop when the root is reached
    if (length(edge) == 0) {
      break
    }

    edge_chain <- c(edge_chain, edge)

    ## Move to the parent node
    parent <- f_tree$edge[edge, 1]

    ## Extant descendants of the parent
    parent_tips <- f_tree$tip.label[
      unlist(Descendants(f_tree, parent, type = "tips"))
    ]

    ## Ignore fossil tips
    parent_tips <- parent_tips[
      !grepl("^f[0-9]+$", parent_tips)
    ]

    ## Stop once we've climbed above the original branch
    if (!setequal(parent_tips, tip_names)) {
      break
    }

    current_node <- parent
  }

  rev(edge_chain)
}

## Divides a stochastic character map into consecutive pieces
## whose lengths match a supplied vector of target branch
## lengths. Used to redistribute one original branch's history
## across the several new edges created when fossils are
## inserted into it.
split_map_by_lengths <- function(map, lengths, tol = 1e-6) {
  if (abs(sum(map) - sum(lengths)) > tol) {
    stop(
      sprintf(
        "Lengths do not sum to map length (map = %.12f, lengths = %.12f)",
        sum(map), sum(lengths)
      )
    )
  }

  result <- vector("list", length(lengths))

  current_map <- map

  for (i in seq_along(lengths)) {
    target <- lengths[i]

    piece <- numeric()

    while (target > tol) {
      seg_length <- current_map[1]
      seg_state <- names(current_map)[1]

      if (seg_length <= target + tol) {
        piece <- c(piece, setNames(seg_length, seg_state))

        target <- target - seg_length

        current_map <- current_map[-1]
      } else {
        piece <- c(piece, setNames(target, seg_state))

        current_map[1] <- seg_length - target

        target <- 0
      }
    }

    result[[i]] <- piece
  }

  result
}


############################################################
## SECTION 2: Core simulation pipeline
##
## These four functions are the actual scientific pipeline:
## simulate a tree -> evolve a trait on it -> generate fossils
## -> reconstruct the "observed" fossil tree.
############################################################

## Simulates a birth-death tree, time-slices it before the
## present (to mimic stopping data collection before the
## "present day" in a real fossil record), and repeats until
## at least `min_extant` extant taxa remain.
simulate_tree <- function(b = 0.18,
                           d = 0.10,
                           n = 55,
                           slice_time = 0.5,
                           min_extant = 50) {
  repeat {
    ## Simulate birth-death tree
    tree <- pbtree(b = b, d = d, n = n)

    ## Slice the tree before the present
    tree <- timeSliceTree(tree, slice_time, plot = FALSE)

    ## Calculate node depths
    node_depths <- node.depth.edgelength(tree)
    tree_height <- max(node_depths)

    ## Count extant taxa
    tip_depths <- head(node_depths, Ntip(tree))
    n_extant <- sum(tip_depths == tree_height)

    ## Accept only sufficiently large trees
    if (n_extant >= min_extant) {
      break
    }
  }

  tree <- ladderize(tree)
  tree$tree_height <- tree_height
  tree$node_depths <- node_depths

  return(tree)
}

## Simulates the complete stochastic history of a binary
## character (states "A"/"B") on a tree, under a continuous-time
## Markov model with a symmetric transition rate `Q`.
simulate_history <- function(tree,
                              Q = matrix(
                                c(-0.05, 0.05,
                                  0.05, -0.05),
                                nrow = 2, byrow = TRUE
                              ),
                              anc = NULL) {
  rownames(Q) <- c("A", "B")
  colnames(Q) <- c("A", "B")

  history <- sim.history(tree, Q = Q, anc = anc)

  tip_states <- getStates(history, type = "tips")
  history$Q <- Q
  return(history)
}

## Simulates a state-dependent fossil record: fossils belonging
## to lineages in state A are preserved at rate `psi_A`, and
## those in state B at rate `psi_B`. This is the core mechanism
## being tested in the study -- when psi_A != psi_B, the fossil
## record is a biased (non-random) sample of the true history.
simulate_fossils <- function(tree,
                              history,
                              slice_time = 0.5,
                              psi_A = 0.2,
                              psi_B = 0.02) {
  ## Tree height and node depths
  node_depths <- node.depth.edgelength(tree)
  tree_height <- max(node_depths)

  ## Simulate fossil occurrences for each character state
  fossils_A <- sim.fossils.poisson(tree = tree, rate = psi_A)
  fossils_B <- sim.fossils.poisson(tree = tree, rate = psi_B)

  fossil_sets <- list(A = fossils_A, B = fossils_B)

  fossils <- data.frame()

  ## Process each preservation regime
  for (state in names(fossil_sets)) {
    fossil_data <- fossil_sets[[state]]

    ## Examine each simulated fossil occurrence
    for (i in seq_along(fossil_data$hmin)) {
      fossil_node <- fossil_data$edge[i]

      ## Convert fossil age to time measured from the root
      fossil_age <- tree_height - fossil_data$hmin[i] + slice_time

      fossil_edge <- which(tree$edge[, 2] == fossil_node)
      parent_node <- tree$edge[fossil_edge, 1]
      branch_start <- node_depths[parent_node]
      distance <- fossil_age - branch_start

      ## True character state at the fossil's exact position
      cumulative_map <- cumsum(history$maps[[fossil_edge]])
      true_state <- names(cumulative_map)[
        match(TRUE, cumulative_map >= distance)
      ]

      ## Keep the fossil only if its preservation "regime"
      ## matches the true state at that point (this is what
      ## makes the fossil record state-dependent rather than
      ## uniform)
      if (true_state == state) {
        fossil <- fossil_data[i, ]
        fossil$state <- state
        fossils <- rbind(fossils, fossil)
      }
    }
  }

  return(fossils)
}

## Inserts each simulated fossil into the tree as a
## zero-length terminal branch, reconstructs the true character
## history on this fossil-inclusive tree, then removes any
## extinct lineages that were never actually sampled as
## fossils. Returns the "observed" fossil tree a real
## palaeontologist would have to work with, plus the TRUE
## state at every tip and node (used later to score accuracy).
build_fossil_tree <- function(tree,
                               history,
                               fossils,
                               slice_time = 0.5) {
  ##########################################################
  ## Insert fossil tips
  ##########################################################
  node_depths <- node.depth.edgelength(tree)
  tree_height <- max(node_depths)

  f_tree <- tree
  fossil_id <- 1

  edits <- data.frame(
    fossil = character(),
    hmin = numeric(),
    state = character(),
    original_child = integer(),
    stringsAsFactors = FALSE
  )
  edits$tips <- list()

  for (child_node in unique(fossils$edge)) {
    branch_fossils <- fossils[fossils$edge == child_node, ]

    ## Insert oldest fossils first
    branch_fossils <- branch_fossils[
      order(branch_fossils$hmin, decreasing = TRUE),
    ]

    descendant_tips <- tree$tip.label[
      unlist(Descendants(tree, child_node, type = "tips"))
    ]

    for (i in seq_len(nrow(branch_fossils))) {
      fossil_age <- tree_height - (branch_fossils$hmin[i] - slice_time)

      if (length(descendant_tips) == 1) {
        current_node <- which(f_tree$tip.label == descendant_tips)
      } else {
        current_node <- mrca.phylo(f_tree, descendant_tips)
      }

      current_height <- node.depth.edgelength(f_tree)[current_node]
      position <- current_height - fossil_age

      row <- nrow(edits) + 1
      edits[row, c("fossil", "hmin", "state", "original_child")] <- list(
        paste0("f", fossil_id),
        branch_fossils$hmin[i],
        branch_fossils$state[i],
        child_node
      )
      edits$tips[[row]] <- descendant_tips

      f_tree <- bind.tip(
        f_tree,
        tip.label = paste0("f", fossil_id),
        edge.length = 0,
        where = current_node,
        position = position
      )

      fossil_id <- fossil_id + 1
    }
  }

  ##########################################################
  ## Reconstruct the stochastic character history
  ##########################################################
  f_history <- history
  f_history$edge <- f_tree$edge
  f_history$edge.length <- f_tree$edge.length
  f_history$tip.label <- f_tree$tip.label
  f_history$Nnode <- f_tree$Nnode
  f_history$maps <- vector("list", nrow(f_tree$edge))

  ## Copy histories for branches unaffected by fossil insertion
  split_children <- unique(edits$original_child)

  for (edge in seq_len(nrow(f_tree$edge))) {
    child <- f_tree$edge[edge, 2]

    ## Skip fossil terminal branches
    if (child <= Ntip(f_tree) && grepl("^f", f_tree$tip.label[child])) {
      next
    }

    descendant_tips <- f_tree$tip.label[
      unlist(Descendants(f_tree, child, type = "tips"))
    ]
    descendant_tips <- descendant_tips[!grepl("^f", descendant_tips)]

    if (length(descendant_tips) == 1) {
      original_child <- which(tree$tip.label == descendant_tips)
    } else {
      original_child <- mrca.phylo(tree, descendant_tips)
    }

    if (!(original_child %in% split_children)) {
      original_edge <- which(tree$edge[, 2] == original_child)
      f_history$maps[[edge]] <- history$maps[[original_edge]]
    }
  }

  ## Split histories on branches containing fossils
  for (child in split_children) {
    branch <- subset(edits, original_child == child)

    original_edge <- which(tree$edge[, 2] == child)
    original_map <- history$maps[[original_edge]]

    edge_chain <- get_edge_chain(f_tree, branch$tips[[1]])
    edge_lengths <- f_tree$edge.length[edge_chain]

    new_maps <- split_map_by_lengths(original_map, edge_lengths)
    f_history$maps[edge_chain] <- new_maps

    ## Fossil tips inherit the state present at the insertion point
    for (i in seq_len(nrow(branch))) {
      fossil_tip <- which(f_tree$tip.label == branch$fossil[i])
      fossil_edge <- which(f_tree$edge[, 2] == fossil_tip)
      f_history$maps[[fossil_edge]] <- setNames(0, branch$state[i])
    }
  }

  ##########################################################
  ## Rebuild simmap attributes (mapped.edge, tip/node states)
  ##########################################################
  stopifnot(!any(sapply(f_history$maps, is.null)))

  states <- sort(unique(names(unlist(f_history$maps))))

  mapped.edge <- matrix(
    0,
    nrow = length(f_history$maps),
    ncol = length(states),
    dimnames = list(NULL, states)
  )

  for (i in seq_along(f_history$maps)) {
    state_lengths <- tapply(
      f_history$maps[[i]],
      names(f_history$maps[[i]]),
      sum
    )
    mapped.edge[i, names(state_lengths)] <- state_lengths
  }
  f_history$mapped.edge <- mapped.edge

  ## Tip states: final state on the branch leading to each tip
  tip_states <- character(Ntip(f_tree))
  for (i in seq_len(Ntip(f_tree))) {
    edge <- which(f_tree$edge[, 2] == i)
    tip_states[i] <- names(tail(f_history$maps[[edge]], 1))
  }
  names(tip_states) <- f_tree$tip.label
  f_history$states <- tip_states

  ## Node states: start/end state of every edge
  edge_start <- vapply(f_history$maps, function(x) names(x)[1], character(1))
  edge_end <- vapply(f_history$maps, function(x) names(x)[length(x)], character(1))
  f_history$node.states <- cbind(edge_start, edge_end)
  storage.mode(f_history$node.states) <- "character"

  ##########################################################
  ## Remove unsampled extinct lineages -> the observed fossil tree
  ##########################################################
  tip_heights <- node.depth.edgelength(f_tree)[1:Ntip(f_tree)]
  tol <- 1e-8

  delete_tips <- which(
    abs(tip_heights - tree_height) > tol &
      !grepl("^f[0-9]+$", f_tree$tip.label)
  )

  f_tree2 <- drop.tip(f_tree, delete_tips)

  ## Recover true tip states on the pruned tree
  tip_states2 <- setNames(character(Ntip(f_tree2)), f_tree2$tip.label)
  for (i in seq_len(Ntip(f_tree2))) {
    old_tip <- which(f_tree$tip.label == f_tree2$tip.label[i])
    old_edge <- which(f_tree$edge[, 2] == old_tip)
    tip_states2[i] <- names(tail(f_history$maps[[old_edge]], 1))
  }

  ## Recover true internal node states on the pruned tree
  node_numbers <- (Ntip(f_tree2) + 1):(Ntip(f_tree2) + f_tree2$Nnode)
  node_states2 <- setNames(character(length(node_numbers)), node_numbers)

  root <- setdiff(f_tree2$edge[, 1], f_tree2$edge[, 2])
  root_edge <- which(f_tree2$edge[, 1] == root)[1]
  root_state <- names(f_history$maps[[root_edge]])[1]

  for (node in node_numbers) {
    if (node == root) {
      node_states2[as.character(node)] <- root_state
      next
    }

    descendant_tips <- f_tree2$tip.label[
      unlist(Descendants(f_tree2, node, type = "tips"))
    ]

    if (length(descendant_tips) == 1) {
      old_child <- which(f_tree$tip.label == descendant_tips)
    } else {
      old_child <- mrca.phylo(f_tree, descendant_tips)
    }

    old_edge <- which(f_tree$edge[, 2] == old_child)
    node_states2[as.character(node)] <- names(tail(f_history$maps[[old_edge]], 1))
  }

  f_tree2$tip_states <- tip_states2
  f_tree2$node_states <- node_states2

  return(f_tree2)
}


############################################################
## SECTION 3: Accuracy metrics
############################################################

## Colless index: a measure of tree (im)balance. Sums the
## absolute difference in tip counts between the two children
## of every internal node -- 0 for a perfectly balanced tree,
## larger for more lopsided trees.
tree_colless <- function(tree) {
  n <- ape::Ntip(tree)
  children <- split(tree$edge[, 2], tree$edge[, 1])

  count_tips <- function(node) {
    if (node <= n) return(1)
    child_nodes <- children[[as.character(node)]]
    sum(vapply(child_nodes, count_tips, numeric(1)))
  }

  internal_nodes <- (n + 1):(n + tree$Nnode)
  colless <- 0

  for (node in internal_nodes) {
    child_nodes <- children[[as.character(node)]]
    if (length(child_nodes) == 2) {
      left_size <- count_tips(child_nodes[1])
      right_size <- count_tips(child_nodes[2])
      colless <- colless + abs(left_size - right_size)
    }
  }

  return(colless)
}

## Internal nodes whose incoming edge has zero length (created
## when a fossil is inserted exactly at a node). Ancestral state
## estimates at these nodes can be degenerate, so they're
## optionally excluded via error_unconstrained below.
zero_length_parent_nodes <- function(tree) {
  if (is.null(tree$edge.length)) stop("Tree has no edge lengths.")
  unique(tree$edge[tree$edge.length == 0, 1])
}

## Mean raw ancestral-state error: for each node, 1 minus the
## predicted probability assigned to the TRUE state, averaged
## across nodes (optionally restricted to a subset of rows).
## This function does the job the rest of the script originally
## called `multRaw()` for -- `multRaw()` was never actually
## defined anywhere (likely typed directly into an R console
## and never saved), but this does the identical calculation
## with `rows = NULL` by default, so it is a safe drop-in
## replacement everywhere `multRaw()` was called.
multRawSubset <- function(prediction, truth, rows = NULL) {
  if (is.matrix(truth)) {
    one_hot <- truth
  } else if (is.atomic(truth)) {
    one_hot <- matrix(0, nrow = length(truth), ncol = ncol(prediction))
    one_hot[cbind(seq_along(truth), truth)] <- 1
  } else {
    stop("truth must be either a matrix or an atomic vector!")
  }

  if (!identical(dim(prediction), dim(one_hot))) {
    stop("prediction and truth have different dimensions.")
  }

  scores <- 1 - rowSums(prediction * one_hot)

  if (!is.null(rows)) scores <- scores[rows]

  mean(scores)
}


############################################################
## SECTION 4: Run one full simulation
##
## Ties SECTION 2 and SECTION 3 together: simulate a tree +
## history + fossil record, reconstruct ancestral states, and
## score the result. Returns a one-row data frame summarising
## everything about that simulation.
############################################################

run_one_simulation <- function(psi_A, psi_B) {
  ## 1. Simulate tree
  tree <- simulate_tree()

  ## 2. Simulate trait history
  history <- simulate_history(tree)

  ## 3. Tree metrics
  tree_height <- max(ape::node.depth.edgelength(tree))
  tree_shape <- tree_colless(tree)

  ## Number of character-state transitions across the whole tree
  state_transitions <- sum(
    vapply(history$maps, function(x) max(0, length(x) - 1), numeric(1))
  )

  ## 4. Simulate fossils with state-dependent preservation
  fossils <- simulate_fossils(tree, history, psi_A = psi_A, psi_B = psi_B)

  ## 5. Build fossil tree
  fossil_tree <- build_fossil_tree(tree, history, fossils)

  ## 6. Fossil and extant tip counts
  fossil_tips <- sum(grepl("^f", fossil_tree$tip.label))
  extant_tips <- sum(!grepl("^f", fossil_tree$tip.label))
  fossil_ratio <- fossil_tips / extant_tips

  ## 7. Ancestral state estimation (equal-rates Mk model)
  ancER <- ace(
    fossil_tree$tip_states,
    fossil_tree,
    type = "discrete",
    model = "ER"
  )

  ## 8-9. Estimated vs true ancestral states
  prediction <- ancER$lik.anc
  true_states <- to.matrix(fossil_tree$node_states, seq = c("A", "B"))

  ## 10. Accuracy score (see SECTION 3 note on multRawSubset)
  error <- multRawSubset(prediction = prediction, truth = true_states)

  ## Error excluding internal nodes with zero-length child edges
  bad_nodes <- zero_length_parent_nodes(fossil_tree)
  node_numbers <- (Ntip(fossil_tree) + 1):(Ntip(fossil_tree) + fossil_tree$Nnode)
  good_rows <- which(!(node_numbers %in% bad_nodes))

  error_unconstrained <- multRawSubset(
    prediction = prediction,
    truth = true_states,
    rows = good_rows
  )

  ## 11. Return results
  return(data.frame(
    psi_A = psi_A,
    psi_B = psi_B,
    fossil_tips = fossil_tips,
    extant_tips = extant_tips,
    fossil_ratio = fossil_ratio,
    error = error,
    tree_shape = tree_shape,
    tree_height = tree_height,
    state_transitions = state_transitions,
    error_unconstrained = error_unconstrained
  ))
}


############################################################
## SECTION 5: Experimental design
##
## 16 scenarios: every combination of 5 heterogeneity levels
## (how different psi_A and psi_B are from each other) crossed
## with 3 overall preservation levels (how large psi_A + psi_B
## is), plus one "no fossils" control where both rates are 0.
############################################################

heterogeneity <- data.frame(
  category = c(
    "None", "Low", "Medium", "High", "Very High",
    "None", "Low", "Medium", "High", "Very High",
    "None", "Low", "Medium", "High", "Very High",
    "No fossils"
  ),
  A = c(
    0.10, 0.08, 0.06, 0.04, 0.02,
    0.25, 0.20, 0.15, 0.10, 0.05,
    0.40, 0.32, 0.24, 0.16, 0.08,
    0.00
  ),
  B = c(
    0.10, 0.12, 0.14, 0.16, 0.18,
    0.25, 0.30, 0.35, 0.40, 0.45,
    0.40, 0.48, 0.56, 0.64, 0.72,
    0.00
  ),
  RateSum = c(
    rep(0.2, 5),
    rep(0.5, 5),
    rep(0.8, 5),
    0
  ),
  OverallPreservation = c(
    rep("Low", 5),
    rep("Medium", 5),
    rep("High", 5),
    "None"
  ),
  stringsAsFactors = FALSE
)


############################################################
## SECTION 6: Run all 8,000 simulations
##
## For each of the 16 scenarios, runs 500 replicate
## simulations. Occasionally a simulation fails for numerical
## reasons (e.g. a degenerate tree) -- these are retried up to
## `max_attempts` times before giving up on that replicate.
## Results are saved to CSV every 100 successful runs so
## nothing is lost if the run is interrupted.
############################################################

res <- data.frame()
max_attempts <- 20
replicates_per_scenario <- 500
output_file <- "simulation_results_revised.csv"

for (h in seq_len(nrow(heterogeneity))) {
  A <- heterogeneity$A[h]
  B <- heterogeneity$B[h]
  category <- heterogeneity$category[h]
  rate_sum <- heterogeneity$RateSum[h]
  overall <- heterogeneity$OverallPreservation[h]

  for (i in seq_len(replicates_per_scenario)) {
    attempt <- 1
    success <- FALSE

    while (!success && attempt <= max_attempts) {
      sim <- tryCatch(
        {
          run_one_simulation(psi_A = A, psi_B = B)
        },
        error = function(e) {
          message(
            sprintf(
              "Simulation failed (%s: A=%.2f, B=%.2f, replicate=%d, attempt=%d): %s",
              category, A, B, i, attempt, e$message
            )
          )
          NULL
        }
      )

      if (!is.null(sim)) {
        success <- TRUE

        sim$category <- category
        sim$A <- A
        sim$B <- B
        sim$RateSum <- rate_sum
        sim$OverallPreservation <- overall
        sim$Replicate <- i

        res <- rbind(res, sim)

        if (nrow(res) %% 100 == 0) {
          write.csv(res, output_file, row.names = FALSE)
          message(paste("Saved", nrow(res), "successful simulations"))
        }
      } else {
        attempt <- attempt + 1
      }
    }

    if (!success) {
      warning(
        sprintf(
          "Giving up after %d attempts (%s: A=%.2f, B=%.2f, replicate=%d)",
          max_attempts, category, A, B, i
        )
      )
    }
  }
}

## Final save
write.csv(res, output_file, row.names = FALSE)
message("All simulations completed and results saved to ", output_file)
