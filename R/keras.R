.default_build_model <- function(input_shape) {
  model <- .do_call_wrap("keras_model_sequential", arglist=list(), pack="keras")
  model <- .do_call_wrap("layer_dense", arglist=list(object=model,units = 64, activation = "relu",
                                                     input_shape = input_shape), pack="keras")
  model <- .do_call_wrap("layer_dense", arglist=list(object=model,units = 64, activation = "relu"), pack="keras")
  model <- .do_call_wrap("layer_dense", arglist=list(object=model,units = 1), pack="keras")
  .do_call_wrap("compile", arglist=list(object=model, optimizer = "rmsprop", loss = "mse", metrics = c("mae")), pack="keras")
}