#' Diagnostic messages when using furrr for parallization
#'
#' @param ...
#' @param domain
#' @param appendLF
#'
#' @return
#' @export
#'
#' @examples
# Define a function called 'immediateMessage' with optional arguments '...'
# 'domain', and 'appendLF', which default to NULL and TRUE respectively
immediateMessage <- function(..., domain = NULL, appendLF = TRUE) {

  # Create a message using the '...' arguments, 'domain', and 'appendLF' options
  msg <- .makeMessage(..., domain = domain, appendLF = appendLF)

  # Get the call that was used to invoke this function
  call <- sys.call()

  # Create a simple message object with the message and call information
  m <- simpleMessage(msg, call)

  # Get the current class(es) of the message object
  cls <- class(m)

  # Remove the "condition" class from the message object (if it exists)
  cls <- setdiff(cls, "condition")

  # Add "immediateCondition" and "condition" classes to the message object
  cls <- c(cls, "immediateCondition", "condition")

  # Update the class(es) of the message object
  class(m) <- cls

  # Print the message to the console
  message(m)

  # Return an invisible version of the message object
  invisible(m)
}
