#' Print figures with dropdown selector
#'
#' Prints HTML code to show a plot based on the value in a dropdown selector.
#' For use with Quarto/Rmarkdown/Knitr, put in a chunk with the option `results='asis'`.
#'
#' @param plot_func A function to produce a plot using a single input from the `selector` input.
#' @param selector One or more values given to `plot_func` to make plots. There will be one plot per input
#' @param label Text used to label the selector dropdown.
#' @param img_class The CSS class used as an ID for the `img` HTML element. This must be unique.
#' @param imglist_class The CSS class used as an ID for the dropdown HTML element. This must be unique.
#'
#' @return
#' @export
#'
#' @examples
#' n <- c(10, 100, 1000)
#' base_plot_func <- function(x) {
#'   hist(rnorm(x))
#' }
#' print_figures_with_selector(base_plot_func, n, 'Count', 'base_test_id')
#'
#' ggplot_func <- function(x) {
#'   df <- data.frame(
#'     sex=factor(rep(c("F", "M"), each=x)),
#'     weight=round(c(rnorm(x, mean=55, sd=5), rnorm(x, mean=65, sd=5)))
#'   )
#'   print(ggplot(df, aes(x=weight)) + geom_histogram())
#' }
#' print_figures_with_selector(ggplot_func, n, 'Number', 'ggplot_test_id')
print_figures_with_selector <- function(plot_func, selector, label, img_class, imglist_class = paste0(img_class, '_list')) {
  # Make plots encoded in base64
  plots <- unlist(lapply(selector, function(x) {
    temp_path <- tempfile(fileext = '.png')
    png(temp_path)
    plot_func(x)
    dev.off()
    output <- base64enc::base64encode(temp_path)
    file.remove(temp_path)
    return(paste0('data:image/png;base64,', output))
  }))

  # Make selector with the base64 plot as value
  cat(paste0('<b>', label, ':  </b>'))
  cat(paste0('<select id="', imglist_class, '">'))
  cat(paste0('  <option value="', plots, '">', n, '</option>', collapse = '\n'))
  cat(paste0('</select>'))

  # Make image showing the first plot
  cat(paste0('<img id="', img_class, '" width="100%" src="', plots[1], '" />'))

  # Add javascript to change with image is shown based on the selector
  cat(paste0('<script type="text/javascript">'))
  cat(paste0('  function setImgSrc(id) {'))
  cat(paste0('    return function (e) {'))
  cat(paste0('      var img = document.getElementById(id);'))
  cat(paste0('      var select = e.target;'))
  cat(paste0('      img.src = select.options[select.selectedIndex].value;'))
  cat(paste0('      return false;'))
  cat(paste0('    };'))
  cat(paste0('  }'))
  cat(paste0('  document.getElementById("', imglist_class, '").onchange = setImgSrc("', img_class, '");'))
  cat(paste0('</script>'))
}

