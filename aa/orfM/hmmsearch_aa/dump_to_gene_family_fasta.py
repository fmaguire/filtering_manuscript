from utils import CARD
if __name__ == '__main__':
    card = CARD('../../../../data/CARD_canonical/card.json')
    card.get_seqs_per_family('protein')
