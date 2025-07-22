.default_build_model <- function(input_shape) {
  keras_model_sequential <- .get_wrap("keras_model_sequential", pack="keras")
  model <- keras_model_sequential()
  layer_dense <- .get_wrap("layer_dense", pack="keras")
  model <- layer_dense(object=model,units = 64, activation = "relu",
                                                     input_shape = input_shape)
  model <- layer_dense(object=model,units = 64, activation = "relu")
  model <- layer_dense(object=model,units = 1)
  compile <- .get_wrap("compile", pack="keras")
  compile(object=model, optimizer = "rmsprop", loss = "mse", metrics = c("mae"))
}