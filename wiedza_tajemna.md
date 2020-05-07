# Przygotowania

Weźmisz czarno kure...

# Podstawowe informacje

Do wysyłania nowych paczek używam narzędzia `twine`.
Można je zaistalować np. za pomocą `pip install twine`.

## Magiczny plik `.pypirc`

Hasła i adresy repozytoriów są trzymane w pliku `.pypirc` w katalogu domowym:

```ini
[distutils]
index-servers=
    pypi
    testpypi

[testpypi]
repository: https://test.pypi.org/legacy/
username: juzer
password: pasłord

[pypi]
#repository: https://pypi.python.org/pypi
username: juzer
password: pasłord
```

Powyższy przykład zawiera dwa repozytoria: normalne i testowe. Do działania wystarczy normalne.
Nazwa użytkownika i hasło powinny odpowiadać temu co jest używane do logowania do [pypi.org](https://pypi.org).

# Rejestracja

Program `twine` można użyć do rejestracji nowego projektu: `twine register aquaduct`. To zostało już dawno zrobione, Aqua-Duct jest już zarejestrowany.

# Budowanie nowej paczki

Jeśli dobrze pamiętam, w pierwszej publikacji zobowiązaliśmy się do udostępniania kodu źródłowego na [pypi.org](https://pypi.org) zawsze buduję paczkę źródłową:

```sh
python setup.py sdist
```

To polecenie spowoduje powstanie nowego archiwum źródłowego.

# Publikowanie paczki

Założenie: `setup.py` jest OK.
Teraz wkracza do akcji twine:

```sh
twine upload ścieżka/do/nowego/archiwum/źródłowego/zbudowanego/przed/chwilą
```

Wułala.

# Dokumentacja

## Budowanie dokumentacji

Teoretycznie sprawa jest prosta. Wystarczy uruchomić dwa skrypty:

1. `make_html.sh`
1. `make_pdf.sh`

W praktyce trzeba to i owo zainstalować: sphinx i przyjaciele, dość kompletny texlive, i może jeszcze parę rzeczy.

## Publikowanie dokumentacji

Jest jeszcze jeden magiczny skrypt: `cp_html_.sh`. Zakłada on, że ścieżka od katalogu `docs` do katalogu z klonem repozytorium TUNNELING-GROUP.github.io to `../../TUNNELING-GROUP.github.io/`.

W tym katalogu jest katalog `aqua-duct` który zawiera dokumentację Aqua-Duct.
Wspomniany magiczny skrypt `cp_html_.sh` magicznie kopiuje nową dokumentację i chyba nawet wysyła ją do repozytorium.

# Oficjalne repozytorium GitHub

Tak jak było to omawiane - może już czas coś zmienić.

Aby jednak wysłać oficjalne zmiany trzeba:

1. Zrobić piękny merdż z gałęzią `v1.0`.
1. Zrobić ładnie odpowiedni tag.
1. Użyć skryptu `fetch_a_new_.sh` który magicznie wysyła wybrane gałęzie i tagi.

# Starsza magia

Niektóre z tych skryptów mogą sprawiać jakieś kłopoty. Adept wiedzy tajemnej powinien najpierw zapoznać się z ich treścią i może wyłączyć polecenia typu `push`.

Przypomniał mi się jeszcze jeden cytat, jakoś tu dziwnie pasuje. Chodzi mi o "Siostrzeńca czarodzieja" C.S. Lewisa:

> Wędrowcze z dalekich stron,
>
> Uderz w dzwon i czekaj na niebezpieczeństwa,
>
> Lub do szaleństwa łam sobie głowę,
>
> Co byś przeżył gdybyś uderzył.

Cytat z głowy więc może nie dokładny ale oddaje ryzyko korzystania z tych skryptów.

