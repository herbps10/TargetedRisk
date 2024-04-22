estimate_riesz_representer_torch <- function(data, baseline, trt, trt_levels, parameter, params) {
  d_out <- 1
  d_in <- length(baseline) + length(trt_levels)
  hidden <- params$hidden
  learning_rate <- params$learning_rate
  epochs <- params$epochs
  dropout <- params$dropout
  torch_seed <- params$seed

  if(!is.na(torch_seed)) torch::torch_manual_seed(torch_seed)

  data_natural <- torch::torch_tensor(cbind(
    as.matrix(data$training[, c(baseline)]),
    model.matrix(~-1 + factor(trt, levels = trt_levels), data = data.frame(trt = data$training[[trt]]))
  ))
  data_natural_valid <- torch::torch_tensor(cbind(
    as.matrix(data$validation[, c(baseline)]),
    model.matrix(~-1 + factor(trt, levels = trt_levels), data = data.frame(trt = data$validation[[trt]]))
  ))

  data_shifted <- lapply(trt_levels, \(trt_level) {
    torch::torch_tensor(cbind(as.matrix(data$training[, baseline]), model.matrix(~-1 + factor(trt, levels = trt_levels), data = data.frame(trt = rep(trt_level, nrow(data$training))))))
  })
  data_shifted_valid <- lapply(trt_levels, \(trt_level) {
    torch::torch_tensor(cbind(as.matrix(data$validation[, baseline]), model.matrix(~-1 + factor(trt, levels = trt_levels), data = data.frame(trt = rep(trt_level, nrow(data$validation))))))
  })

  riesz <- torch::nn_sequential(
    torch::nn_linear(length(baseline), hidden),
    torch::nn_elu(),
    torch::nn_linear(hidden, hidden),
    torch::nn_elu(),
    torch::nn_dropout(dropout)
  )

  heads <- lapply(trt_levels, \(trt_level) {
    torch::nn_sequential(
      torch::nn_linear(ifelse(parameter == "smr", hidden, hidden + 1), hidden),
      torch::nn_elu(),
      torch::nn_linear(hidden, hidden),
      torch::nn_elu(),
      torch::nn_dropout(dropout),
      torch::nn_linear(hidden, d_out),
      torch::nn_softplus()
    )
  })

  #Map(\(x) torch::nn_init_normal_(x, 0, 0.1), riesz$parameters)

  optimizer <- torch::optim_adam(
    params = c(riesz$parameters, unlist(lapply(heads, \(head) head$parameters))),
    lr = learning_rate,
    weight_decay = 0
  )

  scheduler <- torch::lr_one_cycle(optimizer, max_lr = learning_rate, total_steps = epochs)
  for (epoch in 1:epochs) {
    shared <- riesz(data_natural[, 1:length(baseline)])
    natural <- lapply(seq_along(trt_levels), \(index) {
      if(parameter == "direct") {
        heads[[index]](torch::torch_cat(list(shared, data_natural[, length(baseline) + index, drop = FALSE]), 2))
      }
      else {
        heads[[index]](shared)
      }
    })

    if(parameter == "direct") {
      shifted <- lapply(seq_along(trt_levels), \(index) {
        heads[[index]](torch::torch_cat(list(shared, data_shifted[[index]][, length(baseline) + index, drop = FALSE]), 2))
      })
    }

    # Losses
    if(parameter == "direct") {
      losses <- lapply(seq_along(trt_levels), \(index) {
        (natural[[index]]$pow(2) - 2 * shifted[[index]])$mean(dtype = torch::torch_float())
      })
    }
    else {
      losses <- lapply(seq_along(trt_levels), \(index) {
        (natural[[index]]$pow(2) - 2 * natural[[index]] * data_natural[, length(baseline) + index, drop = FALSE])$mean(dtype = torch::torch_float())
      })
    }
    loss <- Reduce(`+`, losses)

    if (epoch %% 20 == 0) {
      cat("Epoch: ", epoch, " Loss: ", loss$item(), "\n")
    }

    optimizer$zero_grad()
    loss$backward()

    optimizer$step()
    scheduler$step()
  }

  if(parameter == "direct") {
    pred <- lapply(seq_along(trt_levels), \(index) {
      heads[[index]](torch::torch_cat(list(riesz(data_shifted_valid[[index]][,1:length(baseline)]), data_shifted_valid[[index]][, length(baseline) + index, drop = FALSE]), 2))
    })
  }
  else {
    shared_valid <- riesz(data_natural_valid[, 1:length(baseline)])
    pred <- lapply(seq_along(trt_levels), \(index) {
      heads[[index]](shared_valid)
    })
  }

  matrix(unlist(lapply(pred, torch::as_array)), ncol = length(trt_levels), byrow = FALSE)
}
