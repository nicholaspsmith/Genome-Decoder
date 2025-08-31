# 23 and Me Genome Decoder
## A simple python program to decode genome files downloaded from 23 and me.

### What is this?
- Did you delete your 23 and Me account but you would like to view the genetic data that you (hopefully) downloaded before you closed your account?
  - This simple Python script will parse the text file they provided you and will produce a reading of your genetic data.

### Instructions
1. Export genome data from 23andme.com
- [Instructions](https://customercare.23andme.com/hc/en-us/articles/212196868-Accessing-Your-Raw-Genetic-Data)
2. Unzip your genome `.zip` file
3. Locate the text file (.txt)

```bash
# macos
python3 decode_genome.py --file /path/to/genome/file.txt

# linux
/usr/bin/python3 --file /path/to/genome/file.txt
```

### Why?
- 23 and Me filed for bankruptcy in early 2025 and many people were worried that their genomic data would be bought by some company with unknown intentions / morals. It seems that this fear has been mitigated for now, as the company was recently sold to a nonprofit called the TTAM Research Institute which was founded by the original owner and purports to maintain the same privacy standards.
- Still, if you already deleted your account (and downloaded your genomic data!) or if you just don't trust big Corporations (non-profit or otherwise), then this utility should make your genome download human-readable.