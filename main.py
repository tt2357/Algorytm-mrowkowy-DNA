import random
import numpy as np
import Levenshtein as lev
import time


def losowe_dna(dlugosc_dna: int) -> str:  # tworzy losowa nic o podanej dlugosci
    dna = ""
    nukleotydy = ["A", "C", "G", "T"]
    for _ in range(dlugosc_dna):
        dna += random.choice(nukleotydy)

    return dna


def stworzSpektrum(nic_DNA: str, dlugosc_oligo: int) -> list[
    str]:  # tworzy liste oligonukleotydow z wczesniej stworzonej losowej nici, nie dodaje duplikatow do spektrum
    temp = ""
    lista_oligonuk = []
    for i in range(len(nic_DNA) - dlugosc_oligo + 1):
        temp = nic_DNA[i: i + dlugosc_oligo]
        if temp not in lista_oligonuk:
            lista_oligonuk.append(temp)
    return lista_oligonuk


def losowy_oligo(
        dlugosc_oligo: int) -> str:  # tworzy losowy oligonukleotyd (uzywane tylko przy dodawaniu bledow pozytywnych)
    oligo = ""
    nukleotydy = ["A", "C", "G", "T"]
    for _ in range(dlugosc_oligo):
        oligo += random.choice(nukleotydy)
    return oligo


def waga(oligo1: str, oligo2: str,
         dlugosc_oligo: int) -> int:  # w tym programie im bardziej oligonukleotydy pokrywaja sie, tym wieksza waga krawedzi grafu
    dlugosc_oligo1 = dlugosc_oligo - 1
    while dlugosc_oligo1 > 0:
        if oligo1[-dlugosc_oligo1:] == oligo2[:dlugosc_oligo1]:
            waga = dlugosc_oligo1
            return waga
        else:
            dlugosc_oligo1 -= 1

    if dlugosc_oligo1 == 0:
        waga = 0
        return waga


def zbuduj_macierz(lista_oligonuk, dlugosc_oligo):
    macierz = np.zeros((len(lista_oligonuk), len(lista_oligonuk)))
    for i in range(len(lista_oligonuk)):
        for j in range(len(lista_oligonuk)):
            if i == j:
                continue
            else:
                macierz[i][j] = waga(lista_oligonuk[i], lista_oligonuk[j], dlugosc_oligo)
    return macierz


def odbuduj_nic(lista_oligonuk: list[str], dlugosc_dna: int,
                macierz):  # obecny oznacza obecny indeks a index oznacza nastepny indeks z grubsza
    odbudowana = ""
    sprawdzanie = []
    sprawdzanie.append(0)
    obecny = 0
    while len(odbudowana) < dlugosc_dna:
        p = np.exp(macierz[obecny]) / np.sum(np.exp(macierz[obecny]))
        index = np.random.choice(len(macierz[obecny]),
                                 p=p)  # wybrany zostaje losowy index, im wieksza waga tym wieksza jest szansa na wybranie jego
        najlepsza_waga = int(macierz[obecny][
                                 index])  # najlepsza_waga oznacza wybrana wage, po prostu wczesniej kod dzialal nieco inaczej i poprawiajac nie zmienialem tej zmiennej bo balem sie ze zmienie dodatkowo cos czego nie powinienem
        if najlepsza_waga <= 1:  # jesli jest mniejsza lub rowna 1 to algorytm cofa sie i losuje jeszcze raz
            continue
        if index in sprawdzanie:
            continue  # jesli juz odwiedzil dany wierzcholek to losuje jeszcze raz
        if obecny == 0:
            odbudowana = lista_oligonuk[0]

        odbudowana = odbudowana + lista_oligonuk[index][najlepsza_waga:]
        obecny = index
        sprawdzanie.append(index)

    return odbudowana


def przejscie_mrowek(lista_oligonuk: list[str], dlugosc_oligo: int, macierz, i_mrowek: int, feromon_szanse,
                     macierz_feromonow, dlugosc_dna):
    jakosci = []
    sciezki = []
    lista_wag = []
    odbudowane = []
    for i in range(i_mrowek):
        jakosc = 0
        sprawdzanie = []
        wagi = [0]
        sprawdzanie.append(0)
        obecny = 0
        odbudowana = ""
        while len(odbudowana) < dlugosc_dna:
            if (random.random() <= feromon_szanse) and (sum(macierz_feromonow[obecny]) != 0):
                p = np.exp(macierz_feromonow[obecny]) / np.sum(np.exp(macierz_feromonow[obecny]))
                index = np.random.choice(len(macierz_feromonow[obecny]),p=p)  # im wiekszy wspolczynnik feromonowy tym wieksza szansa na wybranie
                najlepsza_waga = int(macierz_feromonow[obecny][index])
                rzeczywista_waga = int(macierz[obecny][index])
                if najlepsza_waga == 0:  # jesli jest rowna 0 to algorytm cofa sie i losuje jeszcze raz (wlasciwie to jest niepotrzebne ale na razie zostawiam jako zabezpieczenie)
                    continue
                if index in sprawdzanie:
                    continue  # jesli juz odwiedzil dany wierzcholek to losuje jeszcze raz
                if macierz[obecny][
                    index] == dlugosc_oligo - 1:  # licznik jakosci bierze pod uwage to ile jest najwiekszych mozliwych wag w danej sciezce
                    jakosc += 1


            else:
                p = np.exp(macierz[obecny]) / np.sum(np.exp(macierz[obecny]))
                index = np.random.choice(len(macierz[obecny]), p=p)
                najlepsza_waga = int(macierz[obecny][index])  # im wieksza waga tym wieksza szansa na wybranie
                rzeczywista_waga = int(macierz[obecny][index])
                if najlepsza_waga <= 1:  # jesli jest mniejsza lub rowna 1 to algorytm cofa sie i losuje jeszcze raz
                    continue
                if index in sprawdzanie:
                    continue  # jesli juz odwiedzil dany wierzcholek to losuje jeszcze raz
                if najlepsza_waga == dlugosc_oligo - 1:
                    jakosc += 1

            if obecny == 0:
                odbudowana = lista_oligonuk[0]

            odbudowana = odbudowana + lista_oligonuk[index][rzeczywista_waga:]  # to jest wlasciwie tylko po to zeby algorytm wiedzial kiedy konczyc przejscie danej mrowki

            obecny = index
            if obecny != 0:
                sprawdzanie.append(index)
            wagi.append(rzeczywista_waga)

        odbudowane.append(odbudowana)
        sciezki.append(sprawdzanie)
        jakosci.append(jakosc)
        lista_wag.append(wagi)

    jakosci = sorted(jakosci)
    kolejnosc = np.argsort(jakosci)
    sciezki = [sciezki[i] for i in kolejnosc]
    lista_wag = [lista_wag[i] for i in kolejnosc]
    odbudowane = [odbudowane[i] for i in kolejnosc]
    sciezki = sciezki[int(len(sciezki) * 0.9):]
    lista_wag = lista_wag[int(len(lista_wag) * 0.9):]
    jakosci = jakosci[int((len(jakosci)) * 0.9):]
    odbudowane = odbudowane[int((len(odbudowane)) * 0.9):]

    return sciezki, jakosci, lista_wag, odbudowane


def ACO(lista_oligonuk: list[str], dlugosc_oligo: int, macierz, dlugosc_dna: int, iloscMrowek: int,
        wspolczynnikParowania: float):
    start = time.time()
    czas_wykonywania = 300
    wsp_parowania = wspolczynnikParowania
    feromon_szanse = 0
    obecnie_najlepsza_sciezka = []
    obecnie_najlepsza_jakosc = 0
    obecnie_najlepsza_wagi = 0
    macierz_feromonow = np.zeros((len(lista_oligonuk), len(lista_oligonuk)))

    while time.time() - start < czas_wykonywania:
        top = przejscie_mrowek(lista_oligonuk, dlugosc_oligo, macierz, iloscMrowek, feromon_szanse, macierz_feromonow,
                               dlugosc_dna)
        lista_sciezek = top[0]
        lista_jakosci = top[1]
        lista_wag = top[2]
        lista_odbudowanych = top[3]

        # AKTUALIZACJA NAJLEPSZEJ JAK DOTAD SCIEZKI
        if lista_jakosci[-1] > obecnie_najlepsza_jakosc:
            obecnie_najlepsza_sciezka = lista_sciezek[-1]
            obecnie_najlepsza_jakosc = lista_jakosci[-1]
            obecnie_najlepsza_wagi = lista_wag[-1]
            obecnie_najlepsza_odbudowana = lista_odbudowanych[-1]

        # PAROWANIE FEROMONOW
        for i in range(len(macierz_feromonow)):
            for j in range(len(macierz_feromonow[i])):
                macierz_feromonow[i][j] = macierz_feromonow[i][j] * (1 - wsp_parowania)

        # AKTUALIZACJA MACIERZY FEROMONOWEJ
        for i in range(len(lista_sciezek)):
            for j in range(len(lista_sciezek[i])):
                if j == 0:
                    macierz_feromonow[0][lista_sciezek[i][j]] += (0.1 * i)
                    index = lista_sciezek[i][j]
                else:
                    macierz_feromonow[index][lista_sciezek[i][j]] += (
                            0.1 * i)  # im lepsza jakosciowo sciezka tym bardziej zwiekszaja sie feromony
                    index = lista_sciezek[i][j]
        if feromon_szanse < 1:
            feromon_szanse += 0.02

        # ZBUDOWANIE NICI NA PODSTAWIE NAJLEPSZEJ SCIEZKI Z ACO
    odbudowana = obecnie_najlepsza_odbudowana

    return odbudowana


#################################################################################

# DANE OGOLNE
dlugosc_dna = int(input("PODAJ DLUGOSC NICI: "))
dlugosc_oligo = int(input("PODAJ DLUGOSC OLIGONUKLEOTYDU: "))
p_bledow = int(input("PODAJ PROCENT BLEDOW: "))

# DANE OGOLNE
# dlugosc_dna: int = 800
# dlugosc_oligo: int = 8
# p_bledow: int = 10

# DANE DLA ACO
iloscMrowek: int = 100
wspolczynnikParowania: float = 0.1

dlugosc_spektrum_idealnego = dlugosc_dna - dlugosc_oligo + 1

nic_DNA = losowe_dna(dlugosc_dna)
print()
print("PIERWOTNA NIC DNA: ")
print(nic_DNA)

# UTWORZENIE LISTY OLIGONUKLEOTYDOW
lista_oligonuk = stworzSpektrum(nic_DNA, dlugosc_oligo)

# print()
# print("LISTA OLIGO")
# print(len(lista_oligonuk))

ilosc_bledow = int(dlugosc_spektrum_idealnego * (p_bledow / 100))

# ZAPAMIETANIE PIERWSZEGO OLIGONUKLEOTYDU
pierwszy_oligo = lista_oligonuk[0]

# USUNIECIE LOSOWYCH OLIGO
negatywneCoutner: int = 0
while len(lista_oligonuk) > dlugosc_spektrum_idealnego - ilosc_bledow:
    tmp = random.choice(lista_oligonuk)
    if tmp == pierwszy_oligo:
        continue
    else:
        lista_oligonuk.remove(tmp)
        negatywneCoutner += 1

# print()
# print("ILOSC OLIGO PO USUNIECIU OLIGO:")
# print(len(lista_oligonuk))
# print("ILOSC USUNIETYCH OLIGO")
# print(negatywneCoutner)


# WSTAWIENIE NOWYCH OLIGONUKLEOTYDOW ZA USUNIETE
pozytywneCounter = 0
while pozytywneCounter < ilosc_bledow:
    nowy_oligo = losowy_oligo(dlugosc_oligo)
    if nowy_oligo in lista_oligonuk:
        continue
    else:
        lista_oligonuk.append(nowy_oligo)
        pozytywneCounter += 1

# print()
# print("ILOSC OLIGO PO DODANIU LOSOWYCH:")
# print(len(lista_oligonuk))

# print("ILOSC DODANYCH OLIGO:")
# print(pozytywneCounter)

# WYMIESZANIE LISTY
lista_oligonuk.sort()

# USUNIECIE PIERWSZEGO OLIGO I DODANIE GO NA POCZATEK LISTY
lista_oligonuk.remove(pierwszy_oligo)
lista_oligonuk.insert(0, pierwszy_oligo)

# TWORZENIE MACIERZY
macierz = zbuduj_macierz(lista_oligonuk, dlugosc_oligo)

print()
print("##################################################### - GENERATOR ROZWIAZAN LOSOWYCH")
print("DLUGOSC SEKWENCJI:", dlugosc_dna, "| PROCENT BLEDOW:", p_bledow, "%", "| DLUGOSC OLIGO:", dlugosc_oligo)
print()

# WYSWIETLENIE ODBUDOWANEJ NICI PRZY UZYCIU ALGORYTMU GENERATORA ROZWIAZAN LOSOWYCH
odbudowana_nic = odbuduj_nic(lista_oligonuk, dlugosc_dna, macierz)
print()
print("ODBUDOWANA NIC PRZY UZYCIU ALGORYTMU GENERATORA ROZWIAZAN LOSOWYCH: ")
print(odbudowana_nic)

# WYSWIETLENIE ODLEGLOSCI LEVENSHTEINA DLA GENERAOTRA INSNTANCJI LOSOWYCH
print()
print("ODLEGLOSC LEVENSHTEINA DLA GENERATORA: ")
odleglosc = lev.distance(odbudowana_nic, nic_DNA)
print(odleglosc)

################################################################################################################


print()
print("##################################################### - ACO")
print("DLUGOSC SEKWENCJI:", dlugosc_dna, "| PROCENT BLEDOW:", p_bledow, "%", "| DLUGOSC OLIGO:", dlugosc_oligo,
      "| ILOSC MROWEK:", iloscMrowek, "| WSPOLCZYNNIK PAROWANIA:", wspolczynnikParowania * 100, "%")
print()

# #WYSWIETLENIE ODBUDOWANEJ NICI PRZY UZYCIU ALGORYTMU MROWKOWEGO
odbudowana_mrowkowym = ACO(lista_oligonuk, dlugosc_oligo, macierz, dlugosc_dna, iloscMrowek, wspolczynnikParowania)
print()
print("ODBUDOWANA NIC PRZY UZYCIU ALGORYTMU MROWKOWEGO:")
print(odbudowana_mrowkowym)

# WYSWIETLENIE ODLEGLOSCI LEVENSHTEINA DLA ALGORYTMU MROWKOWEGO
print()
print("ODLEGLOSC LEVENSHTEINA DLA MROWKOWEGO: ")
odleglosc = lev.distance(odbudowana_mrowkowym, nic_DNA)
print(odleglosc)