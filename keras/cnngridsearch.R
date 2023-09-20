PARAMS <- flags(
  flag_numeric("epoch", 100),
  flag_numeric("filter", 32),
  flag_numeric("nodes1", 64),
  flag_numeric("nodes2", 32),
  flag_numeric("stride", 3),
  flag_numeric("kernelsize", 3),
  flag_numeric("poolsize", 2),
  flag_string("optimizer", "rmsprop"),
  flag_numeric("lr_annealing", 0.01),
  flag_string("activation1", "relu"),
  flag_string("activation2", "relu")
)
modelCNN <- keras_model_sequential() %>% 
  layer_conv_1d(filters=PARAMS$filter, kernel_size=PARAMS$kernelsize, 
    strides=PARAMS$stride, input_shape=c(nFeature, 1)) %>%
   layer_max_pooling_1d(pool_size=PARAMS$poolsize) %>%
   layer_flatten() %>%
   layer_dense(units=PARAMS$nodes1, activation=PARAMS$activation1) %>%
   layer_dense(units=PARAMS$nodes2, activation=PARAMS$activation2) %>%
   layer_dense(units=1) %>% 
  compile(
    loss = 'mse',
    metrics = c('mean_absolute_error'),
    optimizer = PARAMS$optimizer
  ) %>%
  fit(
    x = testX,
    y = testY,
    epochs = PARAMS$epoch,
    batch_size = 32,
    validation_split = 0.2,
    callbacks = list(
      callback_early_stopping(patience = 5),
      callback_reduce_lr_on_plateau(factor = PARAMS$lr_annealing)
    ),
    verbose = FALSE
  )