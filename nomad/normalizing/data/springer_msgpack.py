# Copyright 2019 Alvin Noe Ladines, Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Generates and queries a msgpack database of springer-related quantities downloaded from
http://materials.springer.com. The database is stuctured as

space_group_number : normalized_formula : springer_id : entry

The msgpack file can be queried using ArchiveFileDB.

The html parser was taken from a collection of scripts from FHI without further testing.
"""

import requests
import re
import os
from bs4 import BeautifulSoup

from nomad.archive_library.filedb import ArchiveFileDB


DB_NAME = '.springer.msg'

spacesRe = re.compile(r"\s+")

symbRe = re.compile(r"[A-Z][a-z]{0,3}")

numRe = re.compile(r"[0-9.]+")

bracketRe = re.compile(r"\[")

closingBraketRe = re.compile(r"\]")

columnNames = {
    "Normalized_formula": "normalized_formula",
    "Alphabetic Formula:": "alphabetic_formula",
    "Classification by Properties:": "classification",
    "Compound Class(es):": "compound_classes",
    "Dataset ID": "id",
    "Space Group:": "space_group_number",
}


def parseSymbol(formulaStr):
    m = symbRe.match(formulaStr)
    if m:
        return (m.group(), formulaStr[len(m.group()):])
    else:
        return (None, formulaStr)


def parseAmount(formulaStr):
    m = numRe.match(formulaStr)
    if m:
        return (float(m.group()), formulaStr[len(m.group()):])
    else:
        return (1.0, formulaStr)


def parseSimpleEntry(formulaStr):
    sym, rest = parseSymbol(formulaStr)
    if sym is None:
        return (None, formulaStr)
    else:
        am, rest = parseAmount(rest)
        res = {}
        res[sym] = am
        return (res, rest)


def parseComplexEntry(formulaStr, flatten=True):
    res = {}
    m = bracketRe.match(formulaStr)
    if m is None:
        return (None, formulaStr)
    else:
        rest = formulaStr[len(m.group()):]
        while True:
            simE, rest = parseEntry(rest)
            if simE is None: break
            if '#' in simE:
                if 'fragments' in res:
                    res['fragments'].append(simE)
                else:
                    res['fragments'] = [simE]
            else:
                for sym, am in simE.items():
                    if sym in res:
                        res[sym] += am
                    else:
                        res[sym] = am
        m2 = closingBraketRe.match(rest)
        if m2 is None:
            return (None, formulaStr)
        rest = rest[len(m2.group()):]
        am, rest = parseAmount(rest)
        for k, v in res.items():
            res[k] = v * am
        else:
            res['#'] = am
        return (res, rest)


def parseEntry(formulaStr):
    e, rest = parseSimpleEntry(formulaStr)
    if e is not None:
        return (e, rest)
    return parseComplexEntry(formulaStr)


def parseFormula(formulaStr):
    res = {}
    rest = formulaStr
    while len(rest) > 0:
        e, rest = parseEntry(rest)
        if e is None:
            raise Exception("could not parse entry from %r, did parse %s and failed with %r" % (formulaStr, res, rest))
        if '#' in e:
            if 'fragments' in res:
                res['fragments'].append(e)
            else:
                res['fragments'] = [e]
        else:
            for sym, am in e.items():
                if sym in res:
                    res[sym] += am
                else:
                    res[sym] = am
    return res


def normalizeFormula(formulaDict):
    oldTot = sum(formulaDict.values())
    res = {}
    for symb, amount in formulaDict.items():
        res[symb] = int(amount / oldTot * 100.0 + 0.5)
    sortedS = list(res.keys())
    sortedS.sort()
    resStr = ""
    for symb in sortedS:
        resStr += symb
        resStr += str(res[symb])
    return resStr


def parse(htmltext):
    """
    Parser the quantities in columnNames from an html text.
    """
    soup = BeautifulSoup(htmltext, "html.parser")
    results = {}
    for el in soup.findAll(attrs={"class": "data-list__content"}):
        for it in el.findAll(attrs={"class": "data-list__item"}):
            key = it.find(attrs={"class": "data-list__item-key"})
            keyStr = key.string
            value = spacesRe.sub(" ", it.find(attrs={"class": "data-list__item-value"}).get_text())
            if value:
                value = value.strip()
            if keyStr:
                keyStr = keyStr.strip()
                if keyStr in columnNames:
                    keyStr = columnNames[keyStr]
                results[keyStr] = value
    if 'classification' in results:
        results['classification'] = [x.strip() for x in results['classification'].split(",")]
    if 'compound_classes' in results:
        results['compound_classes'] = [x.strip() for x in results['compound_classes'].split(",")]
    if 'alphabetic_formula' in results:
        try:
            f = parseFormula(results['alphabetic_formula'])
            normalized_formula = normalizeFormula(f)
        except Exception:
            normalized_formula = None
    results['normalized_formula'] = normalized_formula
    return results


def download_entries(formula, space_group_number):
    """
    Downloads the springer quantities related to a structure from springer.
    """
    entries = {}
    root = 'https://materials.springer.com/textsearch?searchTerm=%s&datasourceFacet=sm_isp&substanceId=' % formula
    response = requests.get(root)
    if response.status_code != 200:
        return entries
    re_search = re.compile(" href=\"(/isp/[^\"]+)")
    paths = re_search.findall(response.text)
    paths = ['http://materials.springer.com%s' % p for p in paths]
    for path in paths:
        response = requests.get(path)
        if response.status_code != 200:
            continue
        try:
            data = parse(response.text)
        except Exception:
            continue
        space_group_number = data.get('space_group_number', None)
        normalized_formula = data.get('normalized_formula', None)
        id = data.get('id', None)
        if space_group_number is None or normalized_formula is None or id is None:
            continue
        aformula = data.get('alphabetic_formula', None)
        compound = data.get('compound_classes', None)
        classification = data.get('classification', None)
        entry = dict(id=id, aformula=aformula, url=path, compound=compound, classification=classification)
        entries = ArchiveFileDB.merge_dict(entries, {str(space_group_number): {normalized_formula: {id: entry}}})
    return entries


def get_springer_data(normalized_formula, space_group_number):
    """
    Queries a msgpack database for springer-related quantities. Downloads data if not
    found in database and adds it to database.
    """
    entries = {}
    mode = 'w'
    if os.path.isfile(DB_NAME):
        db = ArchiveFileDB(DB_NAME, 'r')
        entries = db.query({str(space_group_number): {normalized_formula: '*'}})
        db.close()
        mode = 'w+'
    if not entries:
        formula = numRe.sub('', normalized_formula)
        entries = download_entries(formula, space_group_number)
        db = ArchiveFileDB(DB_NAME, mode, 3)
        for key, entry in entries.items():
            db.add_data({key: entry})
        db.close()
    db_dict = {}
    entries = entries.get(str(space_group_number), {}).get(normalized_formula, {})
    for id, entry in entries.items():
        db_dict[id] = {
            'spr_id': id,
            'spr_aformula': entry['aformula'],
            'spr_url': entry['url'],
            'spr_compound': entry['compound'],
            'spr_classification': entry['classification']}
    return db_dict
