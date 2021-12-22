# Review of linux commands

Some of the content overlaps with previous lectures by Masa and others.

## Log onto server

Reminder: our server runs on a Linux operating system.  
The following commands are to 'interact' with that operating system and do NOT work on Windows. Most will work exactly the same on a Mac. Masa's instruction is available [here](https://gist.github.com/masaomi/999d1177c00116e61909220c1d40e32e)

```bash
$ ssh your_BFabric_account_name@172.23.30.6
your_BFabric_account_name@172.23.30.6's password:
```

## Unix commands and working with file types common in bioinformatics

By default, the server uses BASH (Bourne Again SHell) to interpret your commands to tell the computer what to do.

### Tips for working quickly in the terminal

- Scrolling at the prompt (where the `$` always is)
  - up arrow goes to previous command
  - down arrow goes to rescent command
  - `Ctrl + a`: Put your cursor at the beginning of the line
  - `Ctrl + e`: Put your cursor at the end of the line
  - tab completion: Press `tab` after you type the first few chalacters of a file or directory name. It should auto-complete if you typed enough characters for it to figure it out.

- Move between directories

  - To go to the parent directory of the current directory ("up" one directory in the directory structure):

    ```bash
    $ cd ..
    ```

  - To go to the most recent directory you were in before the current one  
    (especially useful if you went from /omg/super/long/path to /here/is/another/long/path in one command):

    ```bash
     $ cd -
    ```

- Want multple windows/panes open at the same time? Try using iTerm2 (you have to download from the internet) or try using something called `tmux` (google it!).

* * *

### `man` pages

These are (usually) good resources for the commands we will run. They include the arguments the unix command expects and what other options are available.

```bash
$ man pwd
$ man less
```

Once the man page is open,

- `d`: page down
- `u`: page up
- `q`: exit page and back to command line prompt.

Google the command and you will find numerous examples of how to use the command plus its various options to do the exact thing you want it to do.

* * *

### Tab delimited files

The output of most bioinformatic software is really just a text file with a lot of lines.  
Many files are output as tab delimited files, in which each column contains a specific type of information.  

For files containing certain types of infomation (raw sequencing reads -> FASTQ; mapped reads -> SAM/BAM; etc), there are standardized formats to which everyone adheres.  

Not only understanding how the files should be formatted, but also knowing how to extract the information you want from them via command line is a quick way to organize/check the output.

- `\t`: tab character
- `\n`: new line

* * *

### I/O stream redirection

1. You can manipulate a data stream using multiple commands on one line using the pipe `|`.

   - Pipe `|` takes the standard output (usually what's printed on the screen) of the first command and immediately enters it as input for the following command.
   - Spaces are not required surrounding the pipe, but it makes it easier to read.

   ```bash
   $ command1 file.txt | command2
   ```

   Some series of commands are just not compatible.
   This should not really be the case for the course, but if you start trying to make the longest single command possible, it just won't work!

2. `>`: redirect standard output (most of what's printed to the screen) to a file and save.

   ```bash
   $ command1 input.txt > output.txt
   ```

3. `>>`: Overwrite any file of the same name every time you run the command.

   ```bash
   $ command1 input2.txt >> output.txt
   ```

* * *

### Merge files using `cat`

- `cat`: conCATenate seceral files
- Often used to quickly display the contents of a file
- prints all the contents of a file as standard output

```bash
# merge
$ cat file1.txt file2.txt > merfed_file.txt

# quick display on screen
$ cat file3 .txt
```

* * *

### View large file with `less`

- `less`: view and search only what you want to see on a screen

  ```bash
  $ less file.txt
  ```

  - Arrow keys: scrolling
  - `/patterns`: search "patterns" further down the document
  - `?patterns`: search "patterns" further up the document
  - `n`: find the next instance of the pattern
  - `q`: exit and back to the command line prompt
  - Useful to combine with pipe

    ```bash
    $ command filename.yxy | less
    ```

<br>

* * *

### Pattern grab with `grep`

- Searches for matches to "pattern" in each file
- Prints the entire line which contains "pattern"
- Case sensitive by default

```bash
grep [-civ] "pattern" file(s).txt
```

`-c`: counts the number of lines it appears in, suppresses printing
`-v`: inverse match: returns the lines that do NOT contain the pattern
`-i`: Perform case insensitive match
`'^word'`: ^ searches for lines beginning only with 'word'
`'word$'`: $ searches for lines only ending with 'word'
`'word\b'`: limit search with \\b (ie "words" not found, "sword" found)

Examples

```bash
$ ln -s /scratch/bio373_2021/data/SNPcalling/reference/Ahal.gff
$ grep "CDS" Ahal.gff | less
$ grep -c "CDS" Ahal.gff
$ grep -vc "CDS" Ahal.gff
```

* * *

### Working with compressed files

The commands above work with uncompressed, plain text files.  
Many tools either output compressed (gzipped, bzipped, etc) files, collaborators will send you compressed files or you should compress your files to save disk space.  
When they are compressed, you need to slightly modify your commands to deal with those files. 

```basj
$ zless file.txt.gz
$ zcat file.txt.gz
$ zdiff file.txt.gz
$ zgrep file.txt.gz
```

These are special. Unfortunately, you can't just put a 'z' in front of any command to have it magically work!

Sometimes, you'll have to pipe commands to make it work:

```bash
$ zcat file.txt.gz | cut -f1 -d "_" > newfile.txt
```

* * *

### Select columns with `cut`

The `cut` command allows you to extract information from specific columns. Downside: you need to know the number(s) of the column(s) you want. Counting starts at 1 from the left most column.

```bash
cut -f [-s]  1,4-6 [-d ","] file.txt
```

`-f`: Select fields (columns); Range or comma separated numbers
`-s`: Return only lines which contain one or more delimiter characters
`-d`: Field delimiter. Tab (`\t`) is default.

Examples

```bash
$ cut -f 1,3-5,7 Ahal.gff | less
$ cut -f 1 –s Ahal.gff  | less
```

* * *

### Count stuff with `wc`

Counts and prints number of lines, words, and bytes (all three by default) for each file.

`-l`: counts only lines
`-c`: counts only bytes (normal ASCII characters are 1 byte)
`-m`: counts only characters

```bash
$ wc -l file.txt
$ wc -c file1.txt file2.txt
```

* * *

### Translate or delete with `tr`

Replaces characters (Char1) with other characters (Char2) and prints to screen. You can also do character/string replacement with something like `sed`, but that requires knowledge of regular expressions. I won't cover that here, but it's not difficult to use.

```
tr [-d] Char1 [Char2]
```

\-d delete Char1 (no Char2), in single quotes

Reminder: New line = `\n`, tab = `\t`, and space = ' '

Examples

```
$ echo "Hello World" | tr ' ' '\n'
$ cat Ahal.gff | tr '\t' ' ' | less
$ cat Medtr.fa | tr -d '\n' | less
```

`tr` command replaces the text **pairwisingly**.
If you write the command like below,

```bash
$ less input.txt | tr 'ABC' 'XYZ' > output
```

All the text in "input.txt" will be converted from `A` to `X`, from `B` to `Y`, and from `C` to `Z`.

* * *

### `sort` lines alphanumerically

```bash
sort [-rnu] file
```

`-r`: Sort in reverse
`-n`: Use numerical sort
`-u`: Return only unique lines

Examples

```bash
$ echo "10 1 12 11 100 2" | tr ' ' '\n' | sort
$ echo "10 1 12 11 100 2" | tr ' ' '\n' | sort -n
$ echo "10 1 12 11 100 2" | tr ' ' '\n' | sort -nr
$ echo "12 10 11 12 12 11" | tr ' ' '\n' | sort -u
```

* * *

### Remove duplicate lines with `uniq`

Find unique instances of strings.
Considers only **consecutive** duplicates. Need to sort first.

```
uniq [-c] file
```

`-c`: Count the number of occurrences for each output line

- 1st column is the count
- 2nd column is the unique string

Examples

```
$ echo "12 10 11 12 12 11" | tr ' ' '\n' | uniq
$ echo "12 10 11 12 12 11" | tr ' ' '\n' | sort | uniq
$ echo "12 10 11 12 12 11" | tr ' ' '\n' | sort | uniq -c
```

* * *

### Save space with symlinks

Often you'll want/need to have a file in directories in several locations. Instead of copying the file everytime, you can create a symbolic link (symlink) to the original file location.

```bash
$ ln -s /path/to/file.txt targetName.txt
$ ln -s /scratch/bio373_2021/data/SNPcalling/reference/Ahal.gff /scratch/bio373_2021/YOUR_USERNAME 
```

* * *
