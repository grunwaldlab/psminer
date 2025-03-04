#' Print figures with dropdown selector
#'
#' Prints HTML code to show a plot based on the value in a dropdown selector.
#' For use with Quarto/Rmarkdown/Knitr, put in a chunk with the option
#' `results='asis'`.
#'
#' @param plot_func A function to produce a plot using a single input from the
#'   `selector` input.
#' @param selector A named list of character vectors, that will be used to
#'   generate the input to `plot_func` to make plots. Plots will be made for
#'   every combination of values in the input. The number of character vectors
#'   should equal the number of arguments taken by `plot_func`.
#' @param id_prefix The prefix added to element IDs to distinguish this plot
#'   from others. This must be unique amoung other calls to this functions in a
#'   single HTML file.
#' @param imglist_class The CSS class used as the prefix for the IDs of the
#'   selector dropdown HTML elements. This must be unique amoung other calls to
#'   this functions in a single HTML file.
#' @param hide_single_selector If `TRUE`, don't show selector dropdown if it
#'   contains only one choice.
#' @param zoom If `TRUE`, add JS code to allow the image to be zoomed.
#' @param zoom_slider If `TRUE`, add a slider to the plot to control the zoom.
#' @param ... Passed to [grDevices::png()]
#'
#' @export
#'
#' @examples
#' n <- c(10, 100, 1000)
#' base_plot_func <- function(x) {
#'   hist(rnorm(x))
#' }
#' print_figures_with_selector(base_plot_func, list(Count = n), 'base_test_id')
#'
#' ggplot_func <- function(x) {
#'   df <- data.frame(
#'     sex=factor(rep(c("F", "M"), each=x)),
#'     weight=round(c(rnorm(x, mean=55, sd=5), rnorm(x, mean=65, sd=5)))
#'   )
#'   print(ggplot(df, aes(x=weight)) + geom_histogram())
#' }
#' print_figures_with_selector(ggplot_func, list(Number = n), 'ggplot_test_id')
print_figures_with_selector <- function(plot_func, selector, id_prefix, imglist_class = paste0(id_prefix, '_list'), hide_single_selector = TRUE, zoom = TRUE, zoom_slider = FALSE, ...) {

  # Get combinations of input parameters
  plot_data <- expand.grid(selector, stringsAsFactors = FALSE)

  # Make plots encoded in base64
  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }
  plot_data$base64_plot <- unlist(lapply(seq_len(nrow(plot_data)), function(i) {
    temp_path <- tempfile(fileext = '.png')
    png(temp_path, ...)
    args <- unname(lapply(plot_data, function(column) {
      if (is.list(column)) {
        return(column[[i]])
      } else {
        return(column[i])
      }
    }))
    quiet(do.call(plot_func, args))
    dev.off()
    if (! file.exists(temp_path)) {
      png(temp_path, ...)
      text_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 4, y = 25, size=8, label = "No Plot.") +
        ggplot2::theme_void()
      print(text_plot)
      dev.off()
    }
    output <- base64enc::base64encode(temp_path)
    file.remove(temp_path)
    return(paste0('data:image/png;base64,', output))

  }))
  # Dont show selectors with only one option
  if (hide_single_selector) {
    selector <- selector[lapply(selector, length) > 1]
  }
  plot_data[names(selector)] <- lapply(plot_data[names(selector)], as.character)
  plot_data$plot_id <- apply(plot_data[names(selector)], MARGIN = 1, paste0, collapse = '-')

  cat('\n<!--html_preserve-->\n')

  # Make selectors
  selector_class_id <- paste0(imglist_class, '-', names(selector))
  selector_class_id <- gsub(selector_class_id, pattern = '[^a-zA-Z0-9-]+', replacement = '-')
  names(selector_class_id) <- names(selector)
  for (selector_id in names(selector)) {
    cat(paste0('<b>', selector_id, ':  </b>\n'))
    cat(paste0('<select id="', selector_class_id[selector_id], '">\n'))
    cat(paste0('  <option value="', selector[[selector_id]], '">', selector[[selector_id]], '</option>', collapse = '\n'))
    cat(paste0('\n</select>\n'))
  }

  # Add zooming widget
  # TODO: modify so that it is not relying on URLs to stuff on the internet
  if (zoom) {
    cat(paste0('<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/zoomist@2/zoomist.css" />\n'))
    cat(paste0('<script type="module">\n'))
    cat(paste0('  import Zoomist from "https://cdn.jsdelivr.net/npm/zoomist@2/zoomist.js"\n'))
    cat(paste0('  new Zoomist(".zoomist-container-', id_prefix, '", {\n'))
    cat(paste0('    maxScale: 10,\n'))
    cat(paste0('    bounds: true,\n'))
    cat(paste0('    slider: ', tolower(zoom_slider), ',\n'))
    cat(paste0('    zoomer: true\n'))
    cat(paste0('  })\n'))
    cat(paste0('</script>\n'))
    cat(paste0('<div class="zoomist-container-', id_prefix, '">\n'))
    cat(paste0('  <div class="zoomist-wrapper">\n'))
    cat(paste0('    <div class="zoomist-image">\n'))
  }

  # Make image showing the plot
  img_elem_id = paste0('img-', id_prefix)
  cat(paste0('      <img id="', img_elem_id, '" width="100%" src="', plot_data$base64_plot[1], '" />\n'))

  # Add zooming widget
  if (zoom) {
    cat(paste0('    </div>\n'))
    cat(paste0('  </div>\n'))
    cat(paste0('</div>\n'))
  }

  # Add javascript to change with image is shown based on the selector
  function_name <- paste0(gsub(id_prefix, pattern = '[^a-zA-Z0-9]+', replacement = '_'), '_setImgSrc')
  cat(paste0('<script type="text/javascript">\n'))
  cat(paste0('  function ', function_name, '(event) {\n'))
  cat(paste0('    var plots = {', paste0('"', plot_data$plot_id, '": "', plot_data$base64_plot, '"', collapse = ', '), '};\n'))
  cat(paste0('    var selectorIds = [', paste0('"', selector_class_id, '"', collapse = ", "), '];\n'))
  cat(paste0('    var img = document.getElementById("', img_elem_id, '");\n'))
  cat(paste0('    var selectors = selectorIds.map((id) => document.getElementById(id));\n'))
  cat(paste0('    var plot_id = selectors.map((selector) => selector.options[selector.selectedIndex].value).join("-");\n'))
  cat(paste0('    img.src = plots[plot_id];\n'))
  cat(paste0('    return false;\n'))
  cat(paste0('  }\n'))
  cat(paste0('  document.getElementById("', selector_class_id, '").onchange = ', function_name, ';', collapse = '\n'))
  cat(paste0('\n</script>\n'))

  cat('\n<!--/html_preserve-->\n')

}


