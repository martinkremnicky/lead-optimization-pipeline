import selfies as sf
import random
import functions as fn

ALPHABET = sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols

def bit_flip(selfies_molecule):
    """
    Generate and swap (SELFIES) molecule's element for another from SELFIES alphabet
    """
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list)-1)
    rnd_sample = random.sample(list(ALPHABET),1)[0]
    selfies_molecule_list.pop(rnd_index)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    return fn.validate(''.join(selfies_molecule_list))

def fragment_insertion(selfies_molecule, fragment_size=3, random_size = True):
    """
    Generate and insert SELFIES fragment into a (SELFIES) molecule 

    ### Parameters 
    fragment_size: size of fragment to be generated
        - if ``random_size=True``, becomes upper size limit of fragment generated,
        else (if `=False`) becomes exact fragment size
    """
    if random_size == True:
        fragment_size = random.randint(1,fragment_size)
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list))
    rnd_sample = random.sample(list(ALPHABET), fragment_size) 
    rnd_sample = ''.join(rnd_sample)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    return fn.validate(''.join(selfies_molecule_list))