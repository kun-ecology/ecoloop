#' Save a list with non-dataframe data into csv/txt file
#'
#' @param list_data a list
#' @param data.name  file name for saving the list e.g., "all.csv"
#' @param ... further arguments for write.table
#'
#' @return
#' @export
#'
#' @examples
saveall <- function(list_data,data.name,...)
{
  if (is.list(list_data)==T)
  {
    saveall.name <- data.name
    saveall <- file(saveall.name, open="a")  #creates a file in append mode
    for (i in seq_along(list_data))
    {
      write.table(names(list_data)[i], file=saveall,
                  sep=",", dec=".",quote=FALSE, col.names=FALSE,
                  row.names=FALSE,...)  #writes the name of the list elements ("A", "B", etc)
      write.table(list_data[[i]], file=saveall,
                  sep=",", dec=".", quote=FALSE, col.names=NA,
                  row.names=TRUE,...)  #writes the data.frames
    }

  }
  else
  {
    print("data to be saved must be a list")
  }
  close(saveall)  #close connection to file.csv
}
