# ITIPs
## Identification TE insertions polymorphisms based on a pan-genome and large-scale resequencing data


# Introduction
Identification TE insertions polymorphisms based on a pan-genome and large-scale resequencing data
## Step 1：identification of reference insertions and deletions.

Using each of the 20 *B. rapa* genomes as the reference and aligned each of the left 19 genome sequneces to the reference genome and call reference insertions and deletions using the [smartie-sv](https://github.com/zeeev/smartie-sv) piepeline.

![image](https://github.com/caixu0518/ITIPs/blob/main/png/Step%201.png)

## Step 2: Determination of TE insertions polymorphisms (TIPs).

![image](https://github.com/caixu0518/ITIPs/blob/main/png/Step%202.png)

## Step 3: Genotype TE insertions in a large-scale *B. rapa* population using short reads.
![image](https://github.com/caixu0518/ITIPs/blob/main/png/Step%203.png)

# Installation
One dependent program: [smartie-sv](https://github.com/zeeev/smartie-sv)
The program does not require additional installation.

```
git clone https://github.com/caixu0518/ITIPs.git
cd ITIPs
chmod u+x *pl 
cd scripts
chmod u+x *
```



