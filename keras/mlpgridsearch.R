PARAMS <- flags(
  # Dugumler
  flag_numeric("nodes1", 256),
  flag_numeric("nodes2", 128),
  flag_numeric("nodes3", 64),
  # Birakma oranlari
  flag_numeric("dropout1", 0.4),
  flag_numeric("dropout2", 0.3),
  flag_numeric("dropout3", 0.2),
  # Ogrenme oranlari
  flag_string("optimizer", "rmsprop"),
  flag_numeric("lr_annealing", 0.1),
  # Aktivasyon fonksiyonlari
  flag_string("activation", "relu")
)
modelMLP <- keras_model_sequential() %>%
  layer_dense(units = PARAMS$nodes1, activation = "relu", 
  input_shape = ncol(trainX)) %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = PARAMS$dropout1) %>%
  layer_dense(units = PARAMS$nodes2, activation = PARAMS$activation) %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = PARAMS$dropout2) %>%
  layer_dense(units = PARAMS$nodes3, activation = PARAMS$activation) %>%
  layer_batch_normalization() %>%
  layer_dropout(rate = PARAMS$dropout3) %>%
  layer_dense(units = 1, activation = "linear") %>%
  compile(
    loss = 'mse',
    metrics = c('mean_absolute_error'),
    optimizer = PARAMS$optimizer
  ) %>%
  fit(
    x = testX,
    y = testY,
    epochs = 100,
    batch_size = 32,
    validation_split = 0.2,
    callbacks = list(
      callback_early_stopping(patience = 5),
      callback_reduce_lr_on_plateau(factor = PARAMS$lr_annealing)
    ),
    verbose = FALSE
  )
