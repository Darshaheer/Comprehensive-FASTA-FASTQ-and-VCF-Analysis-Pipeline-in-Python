# -------------------Comprehensive FASTA/FASTQ/VCF Parsing ---------------------------------
# ------------------- A. FASTA Parsing --------------------------------------

with open("sequence.txt", "r") as file1:
    lines1 = file1.readlines()

header = lines1[0].strip()
sequence = ""
for line in lines1[1:]:
    sequence += line.strip()

g_count = sequence.count("G")
c_count = sequence.count("C")
gc_count = ((g_count + c_count)/len(sequence))* 100

Motif = "TATA"
pos_tata = []
for i in range(len(sequence) - len(Motif) + 1):
    part = sequence[i: i + len(Motif)]
    if part == Motif:
        pos_tata.append(i + 1)

print(f"The Header of the Fasta is: {header}")
print(f"The Sequence of the FASTA is: {sequence}")
print(f"The GC content of the Sequence is: {gc_count}")
print(f"The positions of the Motif(TATA) in the sequence is: {pos_tata}")


# ------------------- B. FASTQ Parsing --------------------------------------
ids = []
Sqc = []
QS = []
PS = []
gc_contents = []
high_quality = 0
total_score_sum = 0
total_score_count = 0
max_avg = 0
max_index = 0

with open("SRR390728_1.fastq","rt") as file:
    lines = file.readlines()

num_reads = len(lines)//4

for i in range(0, len(lines), 4):
    ids.append(lines[i].strip())
    sequencee = lines[i + 1].strip()
    Sqc.append(sequencee)
    QS.append(lines[i + 3].strip())

    score = [ord(char) - 33 for char in lines[i + 3].strip()]
    PS.append(score)

    avg_score = sum(score)/len(score)
    if avg_score >= 30:
        high_quality += 1
    total_score_sum += 1
    total_score_count += 1

    g_counts = sequencee.count("G")
    c_counts = sequencee.count("C")
    gc_per = ((g_counts + c_counts)/len(sequencee))* 100
    gc_contents.append(gc_per)

overall_avg = total_score_sum/total_score_count
avg_gc = sum(gc_contents)/len(gc_contents)

from collections import Counter
nucleotide_counts = Counter()
for seq in Sqc:
    nucleotide_counts.update(seq)

for i in range(len(PS)):
    avg = sum(PS[i])/len(PS[i])
    if avg > max_avg:
        max_avg = avg
        max_index = i + 1

print(f"The ID of the first read is: {ids[0]}")
print(f"The Sequence of the first read is: {Sqc[0]}")
print(f"The Quality scores of the first read are: {QS[0]}")
print(f"The phred scores of the first read are: {PS[0]}")
print(f"The overall average phred score is: {overall_avg}")
print(f"Out of total, the high quality reads are: {high_quality}")
print(f"Read number {max_index} ha the highet phred score of {max_avg}")
print(f"The average GC content of the reads is: {avg_gc:.2f}")
print("The GC content of the first 10 reads is:")
for i in range(10):
    print(f"Read number {i + 1} - GC content: {gc_contents[i]}")
print("The reads with low gc content is:")
low_gc = 0
for i in range(len(gc_contents)):
    if gc_contents[i] < 40:
        print(f"Read number {i + 1} - GC content: {gc_contents[i]}")
        low_gc += 1
        if low_gc == 10:
            break

print("The nucleotide requencies are:")
print(f"A: {nucleotide_counts['A']}")
print(f"T: {nucleotide_counts['T']}")
print(f"G: {nucleotide_counts['G']}")
print(f"C: {nucleotide_counts['C']}")

import re
post_hi = []
for i in range(100):
    rec = Sqc[i]
    if re.search(r'TATA',rec):
        post_hi.append(i)
rad = len(post_hi)
print(f"Total number of reads out of first 100, that contin a TATA box is: {rad}")

count = 0
for i in range(len(Sqc)):
    recs = Sqc[i]
    if re.search(r'A{5}',recs):
        print(f"Read number {i + 1} has 5 Adenine's in its sequence: {recs}")
        count += 1
        if count == 10:
            break

for i in range(50):
    nam = Sqc[i]
    if re.search(r'AAA$',nam):
        print(f"Read number {i + 1} has a Sequence that ends with Triple Adenine's: {nam}")

post_tri = []
count2 = 0
for i in range(len(Sqc)):
    usp = Sqc[i]
    if re.search(r'TATA',usp):
        post_tri.append(i + 1)
        count2 += 1
        if count2 == 100:
            break
print(f"The indexes of the reads(first 100) that contain TATA box are: {post_tri}")

cg_total = 0
for i in range(20):
    san = Sqc[i]
    rows = re.findall(r'CG',san)
    cg_total += len(rows)
print(f"The Motif CG has occured {cg_total} times in the first 20 reads")

post_fo = []
for i in range(200):
    chp = Sqc[i]
    ton = re.findall(r'TATA',chp)
    if len(ton) > 1:
        post_fo.append(i)
rede = len(post_fo)
print(f"The number of reads from first 200, that has more than one occurence of TATA box are: {rede}")

count3 = 0
for i in range(len(Sqc)):
    record = Sqc[i]
    if not re.findall(r'G',record):
        print(f"Read number {i + 1} has no Guanine in its sequence: {record}")
        count3 += 1
        if count3 == 5:
            break

post_fi = []
for i in range(len(Sqc)):
    records = Sqc[i]
    if len(Sqc) % 3 == 0 and re.findall(r'TATA',records):
        post_fi.append(i)
reda = len(post_fi)
print(f"Total reads that have a sequence lenght divisible by 3 and have a TATA sequence are: {reda}")

max_len = 0
long_seq = ""
index = 0
for i in range(len(Sqc)):
    leng = Sqc[i]
    if len(leng) > max_len:
        max_len = len(Sqc)
        long_seq = leng
        index = i + 1
print(f"Read number {index} has the longest sequence lenght of {max_len}: {long_seq}")

max_gc = 0
gc_seq = ""
gc_inde = 0
for i in range(len(gc_contents)):
    doc = gc_contents[i]
    docx = re.findall(r'GC',doc)
    if len(soul_king) > max_gc:
        max_gc = len(docx)
        gc_seq = doc
        gc_inde = i + 1
print(f"Read number {gc_inde} has highest gc count of {max_gc}: {gc_seq}")

palindromic_counts = 0
for i in range(500):
    file = Sqc[i]
    if file == file[::-1]:
        palindromic_counts += 1
if palindromic_counts > 1:
    print(f"Total plaindromic sequences found in the first 500 reads are: {palindromic_counts}")
else:
    print("No Palindromic sequences found")


# ------------------- C. VCF Parsing --------------------------------------
count6 = 0
af_count = 0
af2_count = 0
af_total = 0
snp_count = 0
ma_count = 0
ins_count = 0
del_count = 0
other_count = 0
high_qual_var = 0
high_qual_af = 0
gen_ran = []
qual_list = []
common_var = []
rare_var = []
Variant_dicts = []

with open("CH22.vcf","rt") as file2:

    for lineee in file2:
        if lineee.startswith('#'):
            continue
        column = lineee.strip().split('\t')
        chrom = column[0]
        pos = column[1]
        idss = column[2]
        ref = column[3]
        alt = column[4]
        qual = column[5]
        fil = column[6]
        info = column[7]
        format1 = column[8]
        for sample in column[9:12]:
            value = sample.split(":")
            sample_data = dict(zip(format1.split(':'),value))

        if count6 < 5:
            print(f"Variant number {count6}: Chromosome {chrom}, Pos: {pos}, ID: {idss}, Ref.Al: {ref}, Alt.Al: {alt}, Filter: {fil}, Information: {info}")
            print(f"Genotype Information for Var number {count6}: Genotype: {sample_data.get('GT','.\.')}, Genotype Quality: {sample_data.get('GQ','NA')}, and Genotype Depth: {sample_data.get('DP','NA')}")

        count6 += 1

        info_dicts = {}
        data = info.split(";")
        for item in data:
            if "=" in item:
                key, values = item.split("=")
                info_dicts[key] = values
         
        af_strn = info_dicts.get('AF','0')
        af = float(af_strn.split(',')[0])
        ac_strn = info_dicts.get('AC','0')
        ac = float(ac_strn.split(',')[0])
        an = float(info_dicts.get('AN','0'))
        dp = float(info_dicts.get('DP','0')) 

        if af_count < 5:
            print(f"Variant number {af_count}: AF: {af}, AC: {ac}, AN: {an}, DP: {dp}")
        
        af_count += 1
        af_total += af
        
        if af > 0.05:
            high_qual_af += 1
            common_var.append(lineee.strip())
        elif af < 0.01:
            rare_var.append(lineee.strip())

        if af2_count < 5 and af > 0.05:
            print(f"Variant number {af2_count}, Type: Common, Frequency: {af}")
        elif af2_count < 5 and af < 0.01:
            print(f"Variant number {af2_count}, Type: Rare, Frequency: {af}")
        
        af2_count += 1

        try:
            qual_score = float(qual)
            qual_list.append((count6,qual_score))
        except ValueError:
            continue
        
        if qual_score >= 100:
            high_qual_var += 1

        pos_int = int(pos)
        if chrom == "22" and 15000000 <= pos_int <= 20000000:
            gen_ran.append(lineee.strip())

        if "," in alt:
            ma_count += 1
        elif len(ref) == 1 and len(alt) == 1:
            snp_count += 1
        elif len(ref) > len(alt):
            del_count += 1
        elif len(ref) < len(alt):
            ins_count += 1
        else:
            other_count += 1
        
    filtered_variants = common_var + rare_var
    try:
        for fv in filtered_variants:
            column2 = fv.strip().split('\t')
            chrom = column2[0]
            pos = column2[1]
            ref = column2[3]
            alt = column2[4]
            qual = column2[5]
            info = column2[7]
            
            info_dicts = {}
            dataa = info.split(";")
            for item in dataa:
                if "=" in item:
                    key, valuess = item.split("=")
                    info_dicts[key] = valuess
            
            af_strnn = info_dicts.get('AF','0')
            aff = float(af_strnn.split(',')[0])
        Variant_dicts.append({
            "chrom" : chrom,
            "pos" : pos,
            "ref" : ref,
            "alt" : alt,
            "qual" : qual,
            "info" : info,
            "aff" : aff
        })
    except Exception as e:
        print(f"There was an error: {e}")

if af_count > 1:
    avg = af_total/af_count
else:
    avg = 0

import json

with open("Variant Dictionary.json","w") as json_file:
    json.dump(variant_dicts,json_file,indent=4)

if qual_list:
    max_qual = max(qual_list,key=lambda x:x[1])
    min_qual = min(qual_list, key=lambda x:x[1])
    avg_qual = sum(q for _,q in qual_list)/len(qual_list)
else:
    max_qual = min_qual = avg_qual = 0

print(f"The Total SNP's are: {snp_count}")
print(f"The Total Multi-allelic counts are: {ma_count}")
print(f"The Total Insertions are: {ins_count}")
print(f"The Total Deletions are: {del_count}")
print(f"Out of total, the number of high quality variants is: {high_qual_var}")
print(f"Out of total, the number of variants with high allele frequency is: {high_qual_af}")
print(f"Among the variants, the number of common variants is: {len(common_var)}")
print(f"Among the variants, the number of rare variants is: {len(rare_var)}")
print(f"The number of variants with genomic ranges between 15000000 - 20000000 are: {len(gen_ran)}")
print(f"The average Allele frequency is: {avg}")
print(f"Variant number {max_qual[0]} has highest quality score of {max_qual[1]}")
print(f"Variant number {min_qual[0]} has the lowest quality score of {min_qual[1]}")
print(f"The average quality score of the variants is: {avg_qual:.2f}")


# ------------------- D. Exporting Parsed Info --------------------------------------
with open("Parsed Genomic File.txt","w") as out:
    out.write("The Summary of FASTA, FASTQ, and VCF Files\n")
    out.write("================================================\n")
    out.write("The FASTA Extraction\n")
    out.write("================================================\n")
    out.write(f"The Header of the Fasta is: {header}\n")
    out.write(f"The Sequence of the FASTA is: {sequence}\n")
    out.write(f"The GC content of the Sequence is: {gc_count}\n")
    out.write(f"The positions of the Motif(TATA) in the sequence is: {pos_tata}\n")
    out.write("===========================================================================\n")
    out.write("The FASTQ File Extraction\n")
    out.write("===========================================================================\n")
    out.write(f"The ID of the first read is: {ids[0]}\n")
    out.write(f"The Sequence of the first read is: {Sqc[0]} \n ")
    out.write(f"The Quality scores of the first read are: {QS[0]} \n ")
    out.write(f"The phred scores of the first read are: {PS[0]} \n ")
    out.write(f"The overall average phred score is: {overall_avg}\n ")
    out.write(f"Out of total, the high quality reads are: {high_quality}\n ")
    out.write(f"Read number {max_index} ha the highet phred score of {max_avg}\n ")
    out.write(f"The average GC content of the reads is: {avg_gc:.2f}\n ")
    out.write("The GC content of the first 10 reads is: \n ")
    for i in range(10):
        out.write(f"Read number {i + 1} - GC content: {gc_contents[i]} \n ")
    out.write("The reads with low gc content is: \n ")
    low_gc = 0
    for i in range(len(gc_contents)):
        if gc_contents[i] < 40:
            out.write(f"Read number {i + 1} - GC content: {gc_contents[i]} \n ")
            low_gc += 1
            if low_gc == 10:
                break

    out.write("The nucleotide requencies are: \n ")
    out.write(f"A: {nucleotide_counts['A']} \n ")
    out.write(f"T: {nucleotide_counts['T']} \n ")
    out.write(f"G: {nucleotide_counts['G']} \n ")
    out.write(f"C: {nucleotide_counts['C']} \n ")

    import re
    post_hi = []
    for i in range(100):
        rec = Sqc[i]
        if re.search(r'TATA',rec):
            post_hi.append(i)
    rad = len(post_hi)
    out.write(f"Total number of reads out of first 100, that contin a TATA box is: {rad}\n")

    count = 0
    for i in range(len(Sqc)):
        recs = Sqc[i]
        if re.search(r'A{5}',recs):
            out.write(f"Read number {i + 1} has 5 Adenine's in its sequence: {recs}\n")
            count += 1
            if count == 10:
                break

    for i in range(50):
        nam = Sqc[i]
        if re.search(r'AAA$',nam):
            out.write(f"Read number {i + 1} has a Sequence that ends with Triple Adenine's: {nam}\n")

    post_tri = []
    count2 = 0
    for i in range(len(Sqc)):
        usp = Sqc[i]
        if re.search(r'TATA',usp):
            post_tri.append(i + 1)
            count2 += 1
            if count2 == 100:
                break
    out.write(f"The indexes of the reads(first 100) that contain TATA box are: {post_tri}\n")

    cg_total = 0
    for i in range(20):
        san = Sqc[i]
        rows = re.findall(r'CG',san)
        cg_total += len(rows)
    out.write(f"The Motif CG has occured {cg_total} times in the first 20 reads\n")

    post_fo = []
    for i in range(200):
        chp = Sqc[i]
        ton = re.findall(r'TATA',chp)
        if len(ton) > 1:
            post_fo.append(i)
    rede = len(post_fo)
    out.write(f"The number of reads from first 200, that has more than one occurence of TATA box are: {rede}\n")

    count3 = 0
    for i in range(len(Sqc)):
        record = Sqc[i]
        if not re.findall(r'G',record):
            out.write(f"Read number {i + 1} has no Guanine in its sequence: {record}\n")
            count3 += 1
            if count3 == 5:
                break

    post_fi = []
    for i in range(len(Sqc)):
        records = Sqc[i]
        if len(Sqc) % 3 == 0 and re.findall(r'TATA',records):
            post_fi.append(i)
    reda = len(post_fi)
    out.write(f"Total reads that have a sequence lenght divisible by 3 and have a TATA sequence are: {reda}\n")

    max_len = 0
    long_seq = ""
    index = 0
    for i in range(len(Sqc)):
        leng = Sqc[i]
        if len(leng) > max_len:
            max_len = len(Sqc)
            long_seq = leng
            index = i + 1
    out.write(f"Read number {index} has the longest sequence lenght of {max_len}: {long_seq}\n")

    max_gc = 0
    gc_seq = ""
    gc_inde = 0
    for i in range(len(gc_contents)):
        doc = gc_contents[i]
        docx = re.findall(r'GC',doc)
        if len(soul_king) > max_gc:
            max_gc = len(docx)
            gc_seq = doc
            gc_inde = i + 1
    out.write(f"Read number {gc_inde} has highest gc count of {max_gc}: {gc_seq}\n")

    palindromic_counts = 0
    for i in range(500):
        file = Sqc[i]
        if file == file[::-1]:
            palindromic_counts += 1
    if palindromic_counts > 1:
        out.write(f"Total plaindromic sequences found in the first 500 reads are: {palindromic_counts}\n")
    else:
        out.write("No Palindromic sequences found")
    out.write("===========================================================================================\n")
    out.write("The Extraction of VCF File\n")
    out.write("===========================================================================================\n")
    out.write("The Common Variants are:\n")
    out.write("====================================================================================================================\n")
    out.write("CHROMOSOME\tPOSITION\tID\tREF.AL\tALT.AL\tQUALITY\tFILTER\tINFORMATION\tFORMAT\tSAMPLE1\tSAMPLE2\tSAMPLE3\n")
    for cm in common_var:
        column3 = cm.strip().split('\t')
        chrom = column3[0]
        pos = column3[1]
        ids = column3[2]
        ref = column3[3]
        alt = column3[4]
        qual = column3[5]
        fil = column3[6]
        info = column3[7]
        format3 = column3[8]
        dat_st = []
        for dati in column3[9:12]:
            valuees = dati.split(':')
            dati_dat = dict(zip(format3.split(':'),valuees))
            gt = dati_dat.get('GT','.\.')
            gq = dati_dat.get('GQ','NA')
            dp = dati_dat.get('DP','NA')
        dat_st.append(f"GT: {gt}, GQ: {gq}, and DP: {dp}")
        da = "\t".join(dat_st)
    out.write(f"{chrom}\t{pos}\t{ids}\t{ref}\t{alt}\t{qual}\t{fil}\t{info}\t{format3}\t{da}\n")
    out.write("=======================================================================================================================\n")
    out.write("=======================================================================================================================\n")
    out.write("The Rare Variants are:\n")
    for ra in common_var:
        column4 = ra.strip().split('\t')
        chrom = column4[0]
        pos = column4[1]
        ids = column4[2]
        ref = column4[3]
        alt = column4[4]
        qual = column4[5]
        fil = column4[6]
        info = column4[7]
        format4 = column4[8]
        dat_sti = []
        for datii in column4[9:12]:
            valueees = datii.split(':')
            datii_dat = dict(zip(format4.split(':'),valueees))
            gT = datii_dat.get('GT','.\.')
            gQ = datii_dat.get('GQ','NA')
            dP = datii_dat.get('DP','NA')
        dat_st.append(f"GT: {gT}, GQ: {gQ}, and DP: {dP}")
        das = "\t".join(dat_sti)
    out.write(f"{chrom}\t{pos}\t{ids}\t{ref}\t{alt}\t{qual}\t{fil}\t{info}\t{format3}\t{das}\n")
    out.write("=======================================================================================================================\n")
    out.write("=======================================================================================================================\n")
    out.write("The other Relevant Information is:\n")
    out.write(f"The Total SNP's are: {snp_count}\n")
    out.write(f"The Total Multi-allelic counts are: {ma_count}\n")
    out.write(f"The Total Insertions are: {ins_count}\n")
    out.write(f"The Total Deletions are: {del_count}\n")
    out.write(f"Out of total, the number of high quality variants is: {high_qual_var}\n")
    out.write(f"Out of total, the number of variants with high allele frequency is: {high_qual_af}\n")
    out.write(f"Among the variants, the number of common variants is: {len(common_var)}\n")
    out.write(f"Among the variants, the number of rare variants is: {len(rare_var)}\n")
    out.write(f"The number of variants with genomic ranges between 15000000 - 20000000 are: {len(gen_ran)}\n")
    out.write(f"The average Allele frequency is: {avg}\n")
    out.write(f"Variant number {max_qual[0]} has highest quality score of {max_qual[1]}\n")
    out.write(f"Variant number {min_qual[0]} has the lowest quality score of {min_qual[1]}\n")
    out.write(f"The average quality score of the variants is: {avg_qual:.2f}\n")
    out.write("=========================================================================================================================\n")
    out.write("The End of the Extraction Summaries")









