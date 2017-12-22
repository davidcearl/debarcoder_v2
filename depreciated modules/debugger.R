###4x12 layout
size <-  100
bc1_levels <- 4
bc2_levels <- 12
mydf <- data.frame(bc1 = sample(1:bc1_levels, size = size, replace = TRUE), 
           bc2 = sample(1:bc2_levels, size = size, replace = TRUE)
)


ncol <- 6
nrow <- 8

A1 <- c(1,1)
origin <- "A7"
grep("^[[:upper:]]{1}", origin, value = TRUE)
row1<- which(LETTERS == gsub("[[:digit:]]+", "", origin))
col1<- as.numeric(gsub("[[:upper:]]{1}", "", origin))
col1

mydf$wellnumber <- mydf$bc1 + (mydf$bc2-1)*mydf$bc1
row1
row1+nrow-1
LETTERS[row1:(row1+nrow-1)]

matrix(paste0(rep(LETTERS[row1:(row1+nrow-1)], each = ncol), rep(col1:(col1+ncol-1), times = nrow)), ncol = ncol, byrow = TRUE)

bc1_row <- TRUE
layout <- matrix(paste0(rep(LETTERS[row1:(row1+nrow-1)], each = ncol), rep(col1:(col1+ncol-1), times = nrow)), ncol = bc2_levels, byrow = bc1_row)


layoutbc1 <- matrix(rep(1:bc1_levels, each = bc2_levels), nrow = nrow, ncol = ncol, byrow = bc1_row)
layoutbc2 <- matrix(rep(1:bc2_levels, times = bc1_levels), nrow = nrow, ncol = ncol, byrow = bc1_row)
layoutbc2


list(1, 2, 3, 4, 5, 6)[[c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)]]

mapply(function(x, y) layout[x, y], mydf$bc1, mydf$bc2)

###6x8 layout

size <-  10
bc1_levels <- 6
bc2_levels <- 8

mydf <- data.frame(bc1 = sample(1:bc1_levels, size = size, replace = TRUE), 
                   bc2 = sample(1:bc2_levels, size = size, replace = TRUE)
)


ncol <- 6
nrow <- 8

A1 <- c(1,1)
mydf
mydf$wellnumber <- mydf$bc1 + (mydf$bc2-1)*mydf$bc1



bc1_row <- FALSE
layout <- matrix(paste0(rep(LETTERS[row1:(row1+nrow-1)], each = ncol),
                        rep(col1:(col1+ncol-1), times = nrow)
                        ),
                 ncol = bc2_levels,
                 byrow = bc1_row)
layout
layoutbc1 <- matrix(rep(1:bc1_levels, each = bc2_levels),
                    nrow = nrow, ncol = ncol, byrow = bc1_row)
layoutbc2 <- matrix(rep(1:bc2_levels, times = bc1_levels),
                    nrow = nrow, ncol = ncol, byrow = bc1_row)


layout_df <- data.frame(row = rep(row1:(row1+nrow-1), times = ncol),
           col = rep(col1:(col1+ncol-1), each = nrow),
           colchar = rep(LETTERS[col1:(col1+ncol-1)], each = nrow),
           bc1 = as.numeric(layoutbc1),
           bc2 = as.numeric(layoutbc2)
           )


ggplot2::ggplot(layout_df, ggplot2::aes(x = col, y = row)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = as.factor(-bc1)), col = "grey50") + 
  ggplot2::scale_x_continuous(breaks = col1:(col1+ncol-1), labels = LETTERS[col1:(col1+ncol-1)], position = "top", expand = c(0,0)) + 
  ggplot2::scale_y_reverse(breaks = row1:(row1+nrow-1), expand = c(0,0)) + 
  ggplot2::scale_fill_brewer(type = "seq", palette = "Oranges", guide = FALSE) + 
  ggplot2::geom_text(ggplot2::aes(label = bc1)) + 
  ggplot2::theme_classic()


ggplot2::ggplot(layout_df, ggplot2::aes(x = col, y = row)) + 
  ggplot2::geom_tile(ggplot2::aes(fill = as.factor(-bc2)), col = "grey50") + 
  ggplot2::scale_x_continuous(breaks = col1:(col1+ncol-1), labels = LETTERS[col1:(col1+ncol-1)], position = "top", expand = c(0,0)) + 
  ggplot2::scale_y_reverse(breaks = row1:(row1+nrow-1), expand = c(0,0)) + 
  ggplot2::scale_fill_brewer(type = "seq", palette = "Blues", guide = FALSE) + 
  ggplot2::geom_text(ggplot2::aes(label = bc2)) + 
  ggplot2::theme_classic()


mydf$well <- mapply(function(x, y) layout[x, y], mydf$bc1, mydf$bc2)
mydf$well
mydf

expID <- 24584
uploadlist <- list.files(file.path(expID, "debarcoded"))
print(uploadlist)
#print(new.exp$id)
i <- 1
file.path(expID, "debarcoded", uploadlist[2])

cyto_session<- CytobankAPI::authenticate("vanderbilt", auth_token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJhODliMGY3M2U4ZTUzNmFiOGNjMDViZDhjMTZhOGQxZCIsImV4cCI6MTUwOTc0NTkxMSwidXNlcl9pZCI6MTQ3LCJhdWQiOiJjeXRvYmFua19hcGlfdjFfdXNlcnMiLCJpYXQiOjE1MDk3MTcxMTEsImlzcyI6Imh0dHBzOi8vdmFuZGVyYmlsdC5jeXRvYmFuay5vcmcvIiwibmJmIjoxNTA5NzE3MTExLCJzdWIiOiJjeXRvYmFua19hcGlfdjEifQ.3fF0IfCtN-vP4h1gf32eyTgazYxpfv_hE0WncauGrUo")
new.exp<- list()
new.exp$id <- expID
for(i in 6:length(uploadlist)) {
  #if (is.function(updateProgress)) {
   # updateProgress(detail = uploadlist[[i]])
  #}
  # i <- 5
  print(i)
  tryCatch( #this causes it to ignore the 
    CytobankAPI::fcs_files.upload(cyto_session, new.exp$id, file.path(expID, "debarcoded", uploadlist[i]), output = "raw")
    , error = function(cond){cond})  
}
?tryCatch
