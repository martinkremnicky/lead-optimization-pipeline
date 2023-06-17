import selfies as sf

ALPHABET = sf.get_semantic_robust_alphabet() # Gets the alphabet of robust symbols
PLACEHOLDER_VALUE = -2
TRANSLATABLE_ALPHABET = [sym for sym in ALPHABET if sf.decoder(sym) != '']
UNTRANSLATABLE_ALPHABET = list(ALPHABET-set(TRANSLATABLE_ALPHABET))