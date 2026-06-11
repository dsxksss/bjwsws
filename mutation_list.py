#!/data/user/yanx/Dev/workdir/miniconda3/bin/python -B
from Bio import PDB
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import three_to_one
import argparse
from anarci import anarci


def get_seqs_from_pdb(pdb_file):
    # 创建PDB解析器
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", pdb_file)

    chain_seqs = {
        chain.id: seq1("".join(residue.resname for residue in chain))
        for chain in structure.get_chains()
    }
    return chain_seqs


def get_res_from_pdb(pdb_file):
    # 创建PDB解析器
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", pdb_file)

    chain_seqs = {
        chain.id: {
            residue.get_id()[1]: three_to_one(residue.resname) for residue in chain
        }
        for chain in structure.get_chains()
    }

    chain_uid_pos = {
        chain.id: {residue.get_id()[1]: i for i, residue in enumerate(chain)}
        for chain in structure.get_chains()
    }

    return chain_seqs, chain_uid_pos


def fetch_chain_type(seq: str, old_chain: str) -> str:
    try:
        _, details, _ = anarci([("id", seq)], scheme="imgt")
        chain_type = details[0][0]["chain_type"]
        return "H" if chain_type == "H" else "L", 1
    except:
        return old_chain, 0  # if error occur, return the old chain type


#
# input: A10-20,A25,A30,B20-30,B35
#
natural_AA = [
    "L",
    "A",
    "G",
    "V",
    "S",
    "E",
    "R",
    "T",
    "I",
    "D",
    "P",
    "K",
    "Q",
    "N",
    "F",
    "Y",
    "M",
    "H",
    "W",
    "C",
]


def single_mut(args):
    wt_chain_seqs = get_seqs_from_pdb(args.pdb)

    chain_seqs, chain_uid_pos = get_res_from_pdb(args.pdb)
    # print(chain_seqs)
    # print(chain_uid_pos)
    chain_mut_res_dict = {}
    tmp = args.muts.strip().split(",")
    # get mutation residues
    for t in tmp:
        t = t.strip()  # A10-20
        chain = t[0]
        t = t[1:]
        res_pos = []
        if "-" in t:
            t = t.split("-")
            left = int(t[0])
            right = int(t[1])
            for i in range(left, right + 1):
                res_pos.append(i)
        else:
            res_pos.append(int(t))

        if chain in chain_mut_res_dict:
            chain_mut_res_dict[chain].extend(res_pos)
        else:
            chain_mut_res_dict[chain] = res_pos
    # print(chain_res_dict)
    omit_AAs = list(args.omit_AAs)
    final_AAs = [aa for aa in natural_AA if aa not in omit_AAs]
    print(final_AAs)

    # mutation list
    mut_H_chains_output = open("Mut_Hchains.fasta", "w")
    mut_L_chains_output = open("Mut_Lchains.fasta", "w")
    cnt = 0
    mut_file_for_foldx = open("mut_list_for_foldx.txt", "w")
    with open(args.outpath, "w") as fout:
        for chain in chain_mut_res_dict:
            wt_seq = wt_chain_seqs[chain]
            print(wt_seq)
            chain_type, is_ab = fetch_chain_type(wt_seq, chain)
            assert (
                is_ab == 1
            ), f"The mutated chain {chain} is not heavy or light chain of antibody"

            for uid in chain_mut_res_dict[chain]:
                for aa in final_AAs:
                    wt = chain_seqs[chain][uid]
                    if wt == aa:
                        continue
                    mutation = f"{wt}{chain}{uid}{aa}"
                    fout.write(f"{mutation}\n")
                    mut_file_for_foldx.write(f"{mutation};\n")
                    cnt += 1

                    # output mutated sequence
                    pos = chain_uid_pos[chain][uid]
                    mut_seq = wt_seq[0:pos] + aa + wt_seq[pos + 1 :]
                    if chain_type == "H":
                        mut_H_chains_output.write(f">{mutation}")
                        mut_H_chains_output.write(f"{mut_seq}\n")
                    elif chain_type in ["L", "K"]:
                        mut_L_chains_output.write(f">{mutation}")
                        mut_L_chains_output.write(f"{mut_seq}\n")
                    else:
                        print(f"{chain} chain type error.")

                # print(cnt)
    mut_file_for_foldx.close()
    mut_H_chains_output.close()
    mut_L_chains_output.close()


def main():
    parser = argparse.ArgumentParser(
        description="Score sequences based on a given structure."
    )
    parser.add_argument(
        "--pdb",
        type=str,
        help="input filepath",
    )
    parser.add_argument(
        "--muts",
        type=str,
        required=True,
        default=None,
        help="A10-20,A25,A30,B20-30,B35, index of UID",
    )

    parser.add_argument(
        "--outpath",
        type=str,
        default="./single_mut_list.txt",
        help="output filepath for scores of variant sequences",
    )
    parser.add_argument(
        "--omit_AAs",
        type=str,
        default="C",
        help="Specify which amino acids should be omitted in the generated sequence, e.g. 'AC' would omit alanine and cystine.",
    )
    # parser.add_argument(
    #     "--min_muts", type=int, default=1, choices=range(1,101),
    #     help="define the minimum number of mutations to scan."
    # )
    # parser.add_argument(
    #     "--max_muts", type=int, default=1, choices=range(1,101),
    #     help="define the maximum number of mutations to scan."
    # )
    # parser.add_argument(
    #     "--max_variants", type=int, default=300000,
    #     help="define the max variants per site. eg. '300000'"
    # )

    args = parser.parse_args()

    single_mut(args)


if __name__ == "__main__":
    main()
