#knit all Rmd in this folder
if(!require('knitr'))
{
	install.packges('knitr')
	library('knitr')
}

files<-list.files()
files.Rmd<-files[grep('.Rmd$',files)]

for(rmdfile in files.Rmd)
{
	knit2html(rmdfile)
}
