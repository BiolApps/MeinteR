#######################################################################
## Package    : MeinteR
## File       : shapes.R
## Functions  : findShapes (exported)
##
## Updated    : 30-04-2019
##
## Title      : Function for the detection of DNA shape features
#######################################################################


#' Find putative conformational DNA changes
#'
#' Predicts conformational changes of DNA shapes, such as
#' minor groove width (MGW), roll, propeller twist (ProT) and helix twist (HeIT)
#' in the unmethylated and methylated context using methyl-DNAshape algorithm.
#'
#' @param bed.data A data frame containing input bed-formatted data
#' @param offset Number of nucleotides expanded in each direction (default:50, max:200)
#' @param shape.plot A logical flag. If TRUE, the function will display a multi-plot with the 
#' conformational changes. (default: FALSE)
#' @return 1/ p-value of the MGW in the unmethylated/methylated CpG context for each sequence
#' @return 2/ p-value of the HeIT in the unmethylated/methylated CpG context for each sequence
#' @return 3/ p-value of the ProT in the unmethylated/methylated CpG context for each sequence
#' @return 4/ p-value of the Roll in the unmethylated/methylated CpG context for each sequence
#' @export

findShapes <- function(bed.data, offset = 50, shape.plot = FALSE) {
  if (offset > 200 | offset < 0) {
    stop('Valid offset: [1..200]')
  }
  if (sum(is.na(bed.data)) > 0) {
    stop("Input dataset contains ", sum(is.na(bed.data)), " missing values.")
  }
  bed.data$strand[is.na(bed.data$strand)] <- "+"
  if (mean(bed.data$end - bed.data$start) <= 1) {
    adjust.C <- 0
  } else {
    adjust.C <- 1
  }
  if (mean(bed.data$end - bed.data$start) > 3) {
    exit("The width of the methylation bed data should be less than 4nt.")
  }
  bed.data$chr <-
    droplevels(bed.data$chr) # Drop levels not observed in the dataset
  gr <-
    GRanges(
      seqnames = bed.data$chr,
      ranges = IRanges(
        start = bed.data$start - offset + adjust.C,
        width = 2 * offset + 1
      )
    )
  
  fn_methy <-
    system.file("extdata", "MethylSample.fa", package = "DNAshapeR")
  suppressMessages(getFasta(
    gr,
    Hsapiens,
    width = 2 * offset + 1,
    filename = fn_methy
  ))
  fap = fn_methy
  m = vector()
  meth.pos = offset + 1
  for (i in seq_len(length(gr))) {
    x = paste0(">seq", i)
    m = rbind(m, x)
    m = rbind(m, meth.pos)
  }
  fn_methy_pos <-
    system.file("extdata", "MethylSamplePos.fa", package = "DNAshapeR")
  
  write.table(
    m,
    file = fn_methy_pos,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  pred <- getShape(fn_methy)
  pred_methy <-
    getShape(fn_methy,
             methylate = TRUE,
             methylatedPosFile = fn_methy_pos)
  #Statistical evaluation
  wt.mgw <- wt.HelT <- wt.ProT <- wt.Roll <- c()
  for (i in 1:nrow(bed.data)) {
    wt <-
      wilcox.test(x = pred$MGW[i, ],
                  y = pred_methy$MGW[i, ],
                  paired = TRUE)
    wt.mgw <- c(wt.mgw, wt$p.value)
    wt <-
      wilcox.test(x = pred$HelT[i, ],
                  y = pred_methy$HelT[i, ],
                  paired = TRUE)
    wt.HelT <- c(wt.HelT, wt$p.value)
    wt <-
      wilcox.test(x = pred$ProT[i, ],
                  y = pred_methy$ProT[i, ],
                  paired = TRUE)
    wt.ProT <- c(wt.ProT, wt$p.value)
    wt <-
      wilcox.test(x = pred$Roll[i, ],
                  y = pred_methy$Roll[i, ],
                  paired = TRUE)
    wt.Roll <- c(wt.Roll, wt$p.value)
  }
  bed.coor <-
    data.frame(paste(bed.data$chr, bed.data$start, bed.data$end, sep = "_"))
  p.MGW <- cbind(bed.coor, wt.mgw)
  colnames(p.MGW) <- c("seqname", "p.MGW")
  p.HelT <- cbind(bed.coor, wt.HelT)
  colnames(p.HelT) <- c("seqname", "p.HelT")
  p.ProT <- cbind(bed.coor, wt.ProT)
  colnames(p.ProT) <- c("seqname", "p.ProT")
  p.Roll <- cbind(bed.coor, wt.Roll)
  colnames(p.Roll) <- c("seqname", "p.Roll")
  if (shape.plot){
  par(mfrow = c(4, 2))
  plotShape(pred$MGW, main = "Unethylated Minor Groove Width")
  panel.first = grid(equilogs = FALSE)
  plotShape(pred_methy$MGW, main = "Methylated Minor Groove Width", colLine = 'red')
  panel.first = grid(equilogs = FALSE)
  plotShape(pred$HelT, main = "Unethylated Helix Twist")
  panel.first = grid(equilogs = FALSE)
  plotShape(pred_methy$HelT, main = "Methylated Helix Twist",  colLine = 'red')
  panel.first = grid(equilogs = FALSE)
  plotShape(pred$ProT, main = "Unethylated Propeller Twist")
  panel.first = grid(equilogs = FALSE)
  plotShape(
    pred_methy$ProT,
    main = "Methylated Propeller Twist",
    ylim = c(-7, -5.5),
    colLine = 'red'
  )
  panel.first = grid(equilogs = FALSE)
  plotShape(pred$Roll, main = "Unethylated Roll")
  panel.first = grid(equilogs = FALSE)
  plotShape(
    pred_methy$Roll,
    main = "Methylated Roll",
    ylim = c(-2, 4),
    colLine = 'red'
  )
  panel.first = grid(equilogs = FALSE)
  dev.off()
  }
  return(list(p.MGW,p.HelT,p.ProT,p.Roll))
}

