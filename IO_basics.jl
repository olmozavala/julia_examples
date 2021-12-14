
# Read a directory
display(readdir())
display(readdir("/home/olmozavala"))
display(readdir(join=true, sort=true)) # Append path

[x for x in readdir() if isdir(x)]# Read folders 
[x for x in readdir() if isfile(x)]# Read files

# Read folders and files recursively
for (root, dirs, files) in  walkdir(".")
    for dir in dirs
        display(dir)
    end
    for file in files
        display(file)
    end
end

# Fuctions
display(joinpath("/home", "olmozavala"))