import selfies as sf
import random



# smiles
benzene = "c1ccccc1"
benzene_sf = sf.encoder(benzene)  # [C][=C][C][=C][C][=C][Ring1][=Branch1]

alphabet=sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols

#def test(variable): #TODO
#    print([ i for i, a in locals().items() if a == variable][0],": ",variable)

#just a 'fancy' printing function
def line(string=''):
    print("-------"+string+"-------")

#converts a (generated) SELFIES molecule to SMILES and back 'validate' it
def validate(selfies_molecule):
    conversion_smi = sf.decoder(selfies_molecule)
    conversion_sf = sf.encoder(conversion_smi)
    return conversion_smi, conversion_sf

# swap molecule's element for another from SELFIES alphabet
def mutate(selfies_molecule):
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list)-1)
    rnd_sample = random.sample(list(alphabet),1)[0]
    selfies_molecule_list.pop(rnd_index)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    return validate(''.join(selfies_molecule_list))

# insert a generated fragment into the molecule 
def mutate_insert(selfies_molecule, fragment_size):
    selfies_molecule_list = list(sf.split_selfies(selfies_molecule))
    rnd_index = random.randint(0,len(selfies_molecule_list))
    rnd_sample = random.sample(list(alphabet), fragment_size) 
    rnd_sample = ''.join(rnd_sample)
    selfies_molecule_list.insert(rnd_index,rnd_sample)
    return validate(''.join(selfies_molecule_list))

#creates derivatives from a single molecule, outputs SMILES list
def generate_derivatives(selfies_molecule,max_fragment_size,total_size):
    finished_derivates = []#[sf.decoder(selfies_molecule)] #add the original in case of population-based algorithms
    for _ in range((total_size/2)):
        fragment_size = random.randint(1,max_fragment_size)
        finished_derivates.append(mutate(selfies_molecule)[0])
        finished_derivates.append(mutate_insert(selfies_molecule,fragment_size)[0])
    return finished_derivates

for i in (generate_derivatives(benzene_sf,10,100)):
    print(i)


















