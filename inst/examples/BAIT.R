\dontrun{

#Not run because this takes a minute or two:
bamPath=file.path(system.file(package='contiBAIT'), 'extdata')

BAIT(bamPath, pairedEnd = FALSE) #our example data is single-end
}
