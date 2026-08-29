#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from flask import Flask, render_template, request, send_file, abort, redirect, url_for, flash
from markupsafe import escape
from werkzeug.utils import secure_filename
from prody.sequence.motif import motifSearch, saveMotifResults, getUniprot
from prody.database.swissprot import SwissProt
from prody.database.refseq import RefSeq
from prody.database.pdb import ProteinDataBank
from prody.utilities.motifhelpers import validateMotif, getLocalDBs, getGenericDBs


app = Flask(__name__)
app.config.from_pyfile("./config.py")
app.static_folder = app.config.get("STATIC_DIR")


@app.route("/", methods=["GET", "POST"])
def index():
    data_dir = app.config.get("DATA_DIR", "./data")
    upload_dir = os.path.join(data_dir, "local")
    os.makedirs(upload_dir, exist_ok=True)
    if request.method == "POST":
        if "file" not in request.files:
            flash("No file was uploaded.")
            return redirect(request.url)
        file = request.files["file"]
        if file.filename == "":
            flash("Please select a file before upload.")
            return redirect(request.url)
        filename = secure_filename(file.filename)
        file.save(os.path.join(upload_dir, filename))
        flash("File {} uploaded.".format(file))
        return redirect(request.url)
    generic_dbs = getGenericDBs(data_dir)
    local_dbs = getLocalDBs(upload_dir)
    return render_template("index.html", generic_dbs=generic_dbs, local_dbs=local_dbs)


@app.route("/searchMotif")
def searchMotif():
    static_dir = app.config.get("STATIC_DIR", "./static")
    data_dir = app.config.get("DATA_DIR", "./data")
    upload_dir = os.path.join(data_dir, "local")
    os.makedirs(upload_dir, exist_ok=True)
    motif = str(escape(request.args.get("motif", default=None, type=str)))
    if not validateMotif(motif):
        flash('"{}" is not PROSITE motif.'.format(motif))
        return redirect(url_for("index"))
    database = str(escape(request.args.get("database", default=None, type=str)))
    try:
        data = motifSearch(database, motif, data_dir)
    except:
        abort(404)
    if data:
        filename = saveMotifResults(database, motif, data, static_dir)
        flash("Results saved to {}".format(filename))
    else:
        flash("No {} pattern found in {}".format(motif, database))
        return redirect(url_for("index"))
    local_dbs = getLocalDBs(upload_dir)
    generic_dbs = getGenericDBs(data_dir)
    return render_template(
        "search.html",
        generic_dbs=generic_dbs,
        local_dbs=local_dbs,
        data=data,
        filename=filename,
        motif=motif,
    )


@app.route("/result/<string:file>")
def download_result(file):
    filename = escape(file)
    download = "{}/{}".format(app.config.get('STATIC_DIR'), filename)
    try:
        flash("File {} has been downloaded.".format(filename))
        return send_file(download, as_attachment=True)
    except:
        abort(404)


@app.route("/download/<string:database>")
def download_database(database):
    db = escape(database)

    def err(*args):
        abort(404)

    {"sp": SwissProt.checkForUpdates, "rs": RefSeq.checkForUpdates, "pdb": ProteinDataBank.downloadRelease,}.get(
        db, err
    )(app.config.get("DATA_DIR"))
    flash("Database {} successfully downloaded!".format(database))
    return redirect(url_for("index"))


@app.route("/alignment/sequence")
def sequence_alignment():
    abort(404)


@app.route("/alignment/structural")
def structural_alignment():
    abort(404)


@app.errorhandler(404)
def page_not_found(error):
    return render_template("404.html")


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5080, debug=True)
