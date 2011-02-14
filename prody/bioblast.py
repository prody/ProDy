#                  Biopython License Agreement
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.

# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE

"""This module file contains Blast module from the Biopython package with
minor modifications."""

"""
Bio/_py3k.py
==============================================================================
"""

# Copyright 2010 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Python 3 compatibility tools (PRIVATE)."""

import sys

if sys.version_info[0] >= 3:
    #Python 3 code (which will be converted using 2to3 script)

    _bytes_to_string = lambda b: b.decode() # bytes to unicode string
    _string_to_bytes = lambda s: s.encode() # unicode string to bytes

    def _as_unicode(s):
        """Turn byte string or unicode string into a unicode string."""
        if isinstance(s, str):
            return s
        #Assume it is a bytes string
        return s.decode()


    def _as_bytes(s):
        """Turn byte string or unicode string into a bytes string.
        
        The Python 2 version returns a (byte) string.
        """
        if isinstance(s, bytes):
            return s
        #Assume it is a unicode string
        return s.encode()
    
    _as_string = _as_unicode

    def _is_int_or_long(i):
        """Check if the value is an integer.

        Note there are no longs on Python 3.
        """
        return isinstance(i, int)

else:
    #Python 2 code

    _bytes_to_string = lambda b: b # bytes to string, i.e. do nothing
    _string_to_bytes = lambda s: str(s) # str (or unicode) to bytes string

    def _as_unicode(s):
        """Turn a (byte) string or a unicode string into a (byte) string."""
        #Will be changed by 2to3 to "isinstance(s, str)" but doesn't matter:
        if isinstance(s, unicode):
            return s
        return s.decode()
    
    def _as_bytes(s):
        """Turn a (byte) string or a unicode string into a (byte) string."""
        return str(s)
    
    _as_string = _as_bytes

    def _is_int_or_long(i):
        """Check if the value is an integer or long."""
        #If the 2to3 long fixer is enabled (which it is by default), this
        #will be changed to "isinstance(i, int) or isinstance(i, int)"
        #but that doesn't matter.
        return isinstance(i, int) or isinstance(i, long)

"""
Bio/Blast/NCBIWWW.py
==============================================================================
"""

# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Patched by Brad Chapman.
# Chris Wroe added modifications for work in myGrid

"""
This module provides code to work with the WWW version of BLAST
provided by the NCBI.
http://blast.ncbi.nlm.nih.gov/

Functions:
qblast        Do a BLAST search using the QBLAST API.
"""

#import sys
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

def qblast(program, database, sequence,
           auto_format=None,composition_based_statistics=None,
           db_genetic_code=None,endpoints=None,entrez_query='(none)',
           expect=10.0,filter=None,gapcosts=None,genetic_code=None,
           hitlist_size=50,i_thresh=None,layout=None,lcase_mask=None,
           matrix_name=None,nucl_penalty=None,nucl_reward=None,
           other_advanced=None,perc_ident=None,phi_pattern=None,
           query_file=None,query_believe_defline=None,query_from=None,
           query_to=None,searchsp_eff=None,service=None,threshold=None,
           ungapped_alignment=None,word_size=None,
           alignments=500,alignment_view=None,descriptions=500,
           entrez_links_new_window=None,expect_low=None,expect_high=None,
           format_entrez_query=None,format_object=None,format_type='XML',
           ncbi_gi=None,results_file=None,show_overview=None, megablast=None,
           ):
    """Do a BLAST search using the QBLAST server at NCBI.

    Supports all parameters of the qblast API for Put and Get.
    Some useful parameters:
    program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
    database       Which database to search against (e.g. "nr").
    sequence       The sequence to search.
    ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
    descriptions   Number of descriptions to show.  Def 500.
    alignments     Number of alignments to show.  Def 500.
    expect         An expect value cutoff.  Def 10.0.
    matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
    filter         "none" turns off filtering.  Default no filtering
    format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
    entrez_query   Entrez query to limit Blast search
    hitlist_size   Number of hits to return. Default 50
    megablast      TRUE/FALSE whether to use MEga BLAST algorithm (blastn only)

    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    http://www.ncbi.nlm.nih.gov/BLAST/blast_overview.html

    """
    import urllib, urllib2
    import time

    assert program in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    # Format the "Put" command, which sends search requests to qblast.
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
    # Additional parameters are taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node9.html on 8 Oct 2010
    parameters = [
        ('AUTO_FORMAT',auto_format),
        ('COMPOSITION_BASED_STATISTICS',composition_based_statistics),
        ('DATABASE',database),
        ('DB_GENETIC_CODE',db_genetic_code),
        ('ENDPOINTS',endpoints),
        ('ENTREZ_QUERY',entrez_query),
        ('EXPECT',expect),
        ('FILTER',filter),
        ('GAPCOSTS',gapcosts),
        ('GENETIC_CODE',genetic_code),
        ('HITLIST_SIZE',hitlist_size),
        ('I_THRESH',i_thresh),
        ('LAYOUT',layout),
        ('LCASE_MASK',lcase_mask),
        ('MEGABLAST',megablast),
        ('MATRIX_NAME',matrix_name),
        ('NUCL_PENALTY',nucl_penalty),
        ('NUCL_REWARD',nucl_reward),
        ('OTHER_ADVANCED',other_advanced),
        ('PERC_IDENT',perc_ident),
        ('PHI_PATTERN',phi_pattern),
        ('PROGRAM',program),
        #('PSSM',pssm), - It is possible to use PSI-BLAST via this API?
        ('QUERY',sequence),
        ('QUERY_FILE',query_file),
        ('QUERY_BELIEVE_DEFLINE',query_believe_defline),
        ('QUERY_FROM',query_from),
        ('QUERY_TO',query_to),
        #('RESULTS_FILE',...), - Can we use this parameter?
        ('SEARCHSP_EFF',searchsp_eff),
        ('SERVICE',service),
        ('THRESHOLD',threshold),
        ('UNGAPPED_ALIGNMENT',ungapped_alignment),
        ('WORD_SIZE',word_size),
        ('CMD', 'Put'),
        ]
    query = [x for x in parameters if x[1] is not None]
    message = urllib.urlencode(query)

    # Send off the initial query to qblast.
    # Note the NCBI do not currently impose a rate limit here, other
    # than the request not to make say 50 queries at once using multiple
    # threads.
    request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi",
                              message,
                              {"User-Agent":"BiopythonClient"})
    handle = urllib2.urlopen(request)

    # Format the "Get" command, which gets the formatted results from qblast
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007    
    rid, rtoe = _parse_qblast_ref_page(handle)
    parameters = [
        ('ALIGNMENTS',alignments),
        ('ALIGNMENT_VIEW',alignment_view),
        ('DESCRIPTIONS',descriptions),
        ('ENTREZ_LINKS_NEW_WINDOW',entrez_links_new_window),
        ('EXPECT_LOW',expect_low),
        ('EXPECT_HIGH',expect_high),
        ('FORMAT_ENTREZ_QUERY',format_entrez_query),
        ('FORMAT_OBJECT',format_object),
        ('FORMAT_TYPE',format_type),
        ('NCBI_GI',ncbi_gi),
        ('RID',rid),
        ('RESULTS_FILE',results_file),
        ('SERVICE',service),
        ('SHOW_OVERVIEW',show_overview),
        ('CMD', 'Get'),
        ]
    query = [x for x in parameters if x[1] is not None]
    message = urllib.urlencode(query)

    # Poll NCBI until the results are ready.  Use a 3 second wait
    delay = 3.0
    previous = time.time()
    while True:
        current = time.time()
        wait = previous + delay - current
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current

        request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi",
                                  message,
                                  {"User-Agent":"BiopythonClient"})
        handle = urllib2.urlopen(request)
        results = _as_string(handle.read())

        # Can see an "\n\n" page while results are in progress,
        # if so just wait a bit longer...
        if results=="\n\n":
            continue
        # XML results don't have the Status tag when finished
        if results.find("Status=") < 0:
            break
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i+len("Status="):j].strip()
        if status.upper() == "READY":
            break

    return StringIO(results)

def _parse_qblast_ref_page(handle):
    """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE).

    The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is proably
    'Request Time of Execution' and RID would be 'Request Identifier'.
    """
    s = _as_string(handle.read())
    i = s.find("RID =")
    if i == -1:
        rid = None
    else:
        j = s.find("\n", i)
        rid = s[i+len("RID ="):j].strip()

    i = s.find("RTOE =")
    if i == -1:
        rtoe = None
    else:
        j = s.find("\n", i)
        rtoe = s[i+len("RTOE ="):j].strip()

    if not rid and not rtoe:
        #Can we reliably extract the error message from the HTML page?
        #e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        #or    "Message ID#32 Error: Query contains no data: Query
        #       contains no sequence data"
        #
        #This used to occur inside a <div class="error msInf"> entry:
        i = s.find('<div class="error msInf">')
        if i != -1:
            msg = s[i+len('<div class="error msInf">'):].strip()
            msg = msg.split("</div>",1)[0].split("\n",1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        #In spring 2010 the markup was like this:
        i = s.find('<p class="error">')
        if i != -1:
            msg = s[i+len('<p class="error">'):].strip()
            msg = msg.split("</p>",1)[0].split("\n",1)[0].strip()
            if msg:
                raise ValueError("Error message from NCBI: %s" % msg)
        #Generic search based on the way the error messages start:
        i = s.find('Message ID#')
        if i != -1:
            #Break the message at the first HTML tag
            msg = s[i:].split("<",1)[0].split("\n",1)[0].strip()
            raise ValueError("Error message from NCBI: %s" % msg)
        #We didn't recognise the error layout :(
        #print s
        raise ValueError("No RID and no RTOE found in the 'please wait' page, "
                         "there was probably an error in your request but we "
                         "could not extract a helpful error message.")
    elif not rid:
        #Can this happen?
        raise ValueError("No RID found in the 'please wait' page."
                         " (although RTOE = %s)" % repr(rtoe))
    elif not rtoe:
        #Can this happen?
        raise ValueError("No RTOE found in the 'please wait' page."
                         " (although RID = %s)" % repr(rid))

    try:
        return rid, int(rtoe)
    except ValueError:
        raise ValueError("A non-integer RTOE found in " \
                         +"the 'please wait' page, %s" % repr(rtoe))

"""
Bio/Blast/NCBIXML.py
==============================================================================
"""
                         
# Copyright 2000 by Bertrand Frottier .  All rights reserved.
# Revisions 2005-2006 copyright Michiel de Hoon
# Revisions 2006-2009 copyright Peter Cock
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""This module provides code to work with the BLAST XML output
following the DTD available on the NCBI FTP
ftp://ftp.ncbi.nlm.nih.gov/blast/documents/xml/NCBI_BlastOutput.dtd

Classes:
BlastParser         Parses XML output from BLAST (direct use discouraged).
                    This (now) returns a list of Blast records.
                    Historically it returned a single Blast record.
                    You are expected to use this via the parse or read functions.

_XMLParser          Generic SAX parser (private).

Functions:
parse               Incremental parser, this is an iterator that returns
                    Blast records.  It uses the BlastParser internally.
read                Returns a single Blast record. Uses the BlastParser internally.
"""
#from Bio.Blast import Record
import xml.sax
from xml.sax.handler import ContentHandler

class _XMLparser(ContentHandler):
    """Generic SAX Parser

    Just a very basic SAX parser.

    Redefine the methods startElement, characters and endElement.
    """
    def __init__(self, debug=0):
        """Constructor

        debug - integer, amount of debug information to print
        """
        self._tag = []
        self._value = ''
        self._debug = debug
        self._debug_ignore_list = []

    def _secure_name(self, name):
        """Removes 'dangerous' from tag names

        name -- name to be 'secured'
        """
        # Replace '-' with '_' in XML tag names
        return name.replace('-', '_')
    
    def startElement(self, name, attr):
        """Found XML start tag

        No real need of attr, BLAST DTD doesn't use them

        name -- name of the tag

        attr -- tag attributes
        """
        self._tag.append(name)
        
        # Try to call a method (defined in subclasses)
        method = self._secure_name('_start_' + name)

        #Note could use try / except AttributeError
        #BUT I found often triggered by nested errors...
        if hasattr(self, method):
            eval("self.%s()" % method)
            if self._debug > 4:
                print "NCBIXML: Parsed:  " + method
        else:
            # Doesn't exist (yet)
            if method not in self._debug_ignore_list:
                if self._debug > 3:
                    print "NCBIXML: Ignored: " + method
                self._debug_ignore_list.append(method)

        #We don't care about white space in parent tags like Hsp,
        #but that white space doesn't belong to child tags like Hsp_midline
        if self._value.strip():
            raise ValueError("What should we do with %s before the %s tag?" \
                             % (repr(self._value), name))
        self._value = ""

    def characters(self, ch):
        """Found some text

        ch -- characters read
        """
        self._value += ch # You don't ever get the whole string

    def endElement(self, name):
        """Found XML end tag

        name -- tag name
        """
        # DON'T strip any white space, we may need it e.g. the hsp-midline
        
        # Try to call a method (defined in subclasses)
        method = self._secure_name('_end_' + name)
        #Note could use try / except AttributeError
        #BUT I found often triggered by nested errors...
        if hasattr(self, method):
            eval("self.%s()" % method)
            if self._debug > 2:
                print "NCBIXML: Parsed:  " + method, self._value
        else:
            # Doesn't exist (yet)
            if method not in self._debug_ignore_list:
                if self._debug > 1:
                    print "NCBIXML: Ignored: " + method, self._value
                self._debug_ignore_list.append(method)
        
        # Reset character buffer
        self._value = ''
        
class BlastParser(_XMLparser):
    """Parse XML BLAST data into a Record.Blast object

    All XML 'action' methods are private methods and may be:
    _start_TAG      called when the start tag is found
    _end_TAG        called when the end tag is found
    """

    def __init__(self, debug=0):
        """Constructor

        debug - integer, amount of debug information to print
        """
        # Calling superclass method
        _XMLparser.__init__(self, debug)
        
        self._parser = xml.sax.make_parser()
        self._parser.setContentHandler(self)
        
        # To avoid ValueError: unknown url type: NCBI_BlastOutput.dtd
        self._parser.setFeature(xml.sax.handler.feature_validation, 0)
        self._parser.setFeature(xml.sax.handler.feature_namespaces, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_pes, 0)
        self._parser.setFeature(xml.sax.handler.feature_external_ges, 0)

        self.reset()

    def reset(self):
        """Reset all the data allowing reuse of the BlastParser() object"""
        self._records = []
        self._header = Header()
        self._parameters = Parameters()
        self._parameters.filter = None #Maybe I should update the class?

    def _start_Iteration(self):
        self._blast = Blast()
        pass

    def _end_Iteration(self):
        # We stored a lot of generic "top level" information
        # in self._header (an object of type Record.Header)
        self._blast.reference = self._header.reference
        self._blast.date = self._header.date
        self._blast.version = self._header.version
        self._blast.database = self._header.database
        self._blast.application = self._header.application

        # These are required for "old" pre 2.2.14 files
        # where only <BlastOutput_query-ID>, <BlastOutput_query-def>
        # and <BlastOutput_query-len> were used.  Now they
        # are suplemented/replaced by <Iteration_query-ID>,
        # <Iteration_query-def> and <Iteration_query-len>
        if not hasattr(self._blast, "query") \
        or not self._blast.query:
            self._blast.query = self._header.query
        if not hasattr(self._blast, "query_id") \
        or not self._blast.query_id:
            self._blast.query_id = self._header.query_id
        if not hasattr(self._blast, "query_letters") \
        or not self._blast.query_letters:
            self._blast.query_letters = self._header.query_letters

        # Hack to record the query length as both the query_letters and
        # query_length properties (as in the plain text parser, see
        # Bug 2176 comment 12):
        self._blast.query_length = self._blast.query_letters
        # Perhaps in the long term we should deprecate one, but I would
        # prefer to drop query_letters - so we need a transition period
        # with both.

        # Hack to record the claimed database size as database_length
        # (as well as in num_letters_in_database, see Bug 2176 comment 13):
        self._blast.database_length = self._blast.num_letters_in_database
        # TODO? Deprecate database_letters next?

        # Hack to record the claimed database sequence count as database_sequences
        self._blast.database_sequences = self._blast.num_sequences_in_database

        # Apply the "top level" parameter information
        self._blast.matrix = self._parameters.matrix
        self._blast.num_seqs_better_e = self._parameters.num_seqs_better_e
        self._blast.gap_penalties = self._parameters.gap_penalties
        self._blast.filter = self._parameters.filter
        self._blast.expect = self._parameters.expect
        self._blast.sc_match = self._parameters.sc_match
        self._blast.sc_mismatch = self._parameters.sc_mismatch

        #Add to the list
        self._records.append(self._blast)
        #Clear the object (a new empty one is create in _start_Iteration)
        self._blast = None

        if self._debug : "NCBIXML: Added Blast record to results"

    # Header
    def _end_BlastOutput_program(self):
        """BLAST program, e.g., blastp, blastn, etc.

        Save this to put on each blast record object
        """
        self._header.application = self._value.upper()

    def _end_BlastOutput_version(self):
        """version number and date of the BLAST engine.

        e.g. "BLASTX 2.2.12 [Aug-07-2005]" but there can also be
        variants like "BLASTP 2.2.18+" without the date.

        Save this to put on each blast record object
        """
        parts = self._value.split()
        #TODO - Check the first word starts with BLAST?

        #The version is the second word (field one)
        self._header.version = parts[1]
        
        #Check there is a third word (the date)
        if len(parts) >= 3:
            if parts[2][0] == "[" and parts[2][-1] == "]":
                self._header.date = parts[2][1:-1]
            else:
                #Assume this is still a date, but without the
                #square brackets
                self._header.date = parts[2]

    def _end_BlastOutput_reference(self):
        """a reference to the article describing the algorithm

        Save this to put on each blast record object
        """
        self._header.reference = self._value

    def _end_BlastOutput_db(self):
        """the database(s) searched

        Save this to put on each blast record object
        """
        self._header.database = self._value

    def _end_BlastOutput_query_ID(self):
        """the identifier of the query

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-ID> is enough
        """
        self._header.query_id = self._value

    def _end_BlastOutput_query_def(self):
        """the definition line of the query

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-def> is enough
        """
        self._header.query = self._value

    def _end_BlastOutput_query_len(self):
        """the length of the query

        Important in old pre 2.2.14 BLAST, for recent versions
        <Iteration_query-len> is enough
        """
        self._header.query_letters = int(self._value)

    def _end_Iteration_query_ID(self):
        """the identifier of the query
        """
        self._blast.query_id = self._value

    def _end_Iteration_query_def(self):
        """the definition line of the query
        """
        self._blast.query = self._value

    def _end_Iteration_query_len(self):
        """the length of the query
        """
        self._blast.query_letters = int(self._value)

##     def _end_BlastOutput_query_seq(self):
##         """the query sequence
##         """
##         pass # XXX Missing in Record.Blast ?

##     def _end_BlastOutput_iter_num(self):
##         """the psi-blast iteration number
##         """
##         pass # XXX TODO PSI

    def _end_BlastOutput_hits(self):
        """hits to the database sequences, one for every sequence
        """
        self._blast.num_hits = int(self._value)

##     def _end_BlastOutput_message(self):
##         """error messages
##         """
##         pass # XXX What to do ?

    # Parameters
    def _end_Parameters_matrix(self):
        """matrix used (-M)
        """
        self._parameters.matrix = self._value
        
    def _end_Parameters_expect(self):
        """expect values cutoff (-e)
        """
        # NOTE: In old text output there was a line:
        # Number of sequences better than 1.0e-004: 1
        # As far as I can see, parameters.num_seqs_better_e
        # would take the value of 1, and the expectation
        # value was not recorded.
        #
        # Anyway we should NOT record this against num_seqs_better_e
        self._parameters.expect = self._value

##     def _end_Parameters_include(self):
##         """inclusion threshold for a psi-blast iteration (-h)
##         """
##         pass # XXX TODO PSI
    
    def _end_Parameters_sc_match(self):
        """match score for nucleotide-nucleotide comparaison (-r)
        """
        self._parameters.sc_match = int(self._value)

    def _end_Parameters_sc_mismatch(self):
        """mismatch penalty for nucleotide-nucleotide comparaison (-r)
        """
        self._parameters.sc_mismatch = int(self._value)

    def _end_Parameters_gap_open(self):
        """gap existence cost (-G)
        """
        self._parameters.gap_penalties = int(self._value)

    def _end_Parameters_gap_extend(self):
        """gap extension cose (-E)
        """
        self._parameters.gap_penalties = (self._parameters.gap_penalties,
                                         int(self._value))

    def _end_Parameters_filter(self):
        """filtering options (-F)
        """
        self._parameters.filter = self._value

##     def _end_Parameters_pattern(self):
##         """pattern used for phi-blast search
##         """
##         pass # XXX TODO PSI

##     def _end_Parameters_entrez_query(self):
##         """entrez query used to limit search
##         """
##         pass # XXX TODO PSI

    # Hits
    def _start_Hit(self):
        self._blast.alignments.append(Alignment())
        self._blast.descriptions.append(Description())
        self._blast.multiple_alignment = []
        self._hit = self._blast.alignments[-1]
        self._descr = self._blast.descriptions[-1]
        self._descr.num_alignments = 0

    def _end_Hit(self):
        #Cleanup
        self._blast.multiple_alignment = None
        self._hit = None
        self._descr = None

    def _end_Hit_id(self):
        """identifier of the database sequence
        """
        self._hit.hit_id = self._value
        self._hit.title = self._value + ' '

    def _end_Hit_def(self):
        """definition line of the database sequence
        """
        self._hit.hit_def = self._value
        self._hit.title += self._value
        self._descr.title = self._hit.title

    def _end_Hit_accession(self):
        """accession of the database sequence
        """
        self._hit.accession = self._value
        self._descr.accession = self._value

    def _end_Hit_len(self):
        self._hit.length = int(self._value)

    # HSPs
    def _start_Hsp(self):
        #Note that self._start_Hit() should have been called
        #to setup things like self._blast.multiple_alignment
        self._hit.hsps.append(HSP())
        self._hsp = self._hit.hsps[-1]
        self._descr.num_alignments += 1
        self._blast.multiple_alignment.append(MultipleAlignment())
        self._mult_al = self._blast.multiple_alignment[-1]

    # Hsp_num is useless
    def _end_Hsp_score(self):
        """raw score of HSP
        """
        self._hsp.score = float(self._value)
        if self._descr.score == None:
            self._descr.score = float(self._value)

    def _end_Hsp_bit_score(self):
        """bit score of HSP
        """
        self._hsp.bits = float(self._value)
        if self._descr.bits == None:
            self._descr.bits = float(self._value)

    def _end_Hsp_evalue(self):
        """expect value value of the HSP
        """
        self._hsp.expect = float(self._value)
        if self._descr.e == None:
            self._descr.e = float(self._value)

    def _end_Hsp_query_from(self):
        """offset of query at the start of the alignment (one-offset)
        """
        self._hsp.query_start = int(self._value)

    def _end_Hsp_query_to(self):
        """offset of query at the end of the alignment (one-offset)
        """
        self._hsp.query_end = int(self._value)

    def _end_Hsp_hit_from(self):
        """offset of the database at the start of the alignment (one-offset)
        """
        self._hsp.sbjct_start = int(self._value)

    def _end_Hsp_hit_to(self):
        """offset of the database at the end of the alignment (one-offset)
        """
        self._hsp.sbjct_end = int(self._value)

##     def _end_Hsp_pattern_from(self):
##         """start of phi-blast pattern on the query (one-offset)
##         """
##         pass # XXX TODO PSI

##     def _end_Hsp_pattern_to(self):
##         """end of phi-blast pattern on the query (one-offset)
##         """
##         pass # XXX TODO PSI

    def _end_Hsp_query_frame(self):
        """frame of the query if applicable
        """
        self._hsp.frame = (int(self._value),)

    def _end_Hsp_hit_frame(self):
        """frame of the database sequence if applicable
        """
        self._hsp.frame += (int(self._value),)

    def _end_Hsp_identity(self):
        """number of identities in the alignment
        """
        self._hsp.identities = int(self._value)

    def _end_Hsp_positive(self):
        """number of positive (conservative) substitutions in the alignment
        """
        self._hsp.positives = int(self._value)

    def _end_Hsp_gaps(self):
        """number of gaps in the alignment
        """
        self._hsp.gaps = int(self._value)

    def _end_Hsp_align_len(self):
        """length of the alignment
        """
        self._hsp.align_length = int(self._value)

##     def _en_Hsp_density(self):
##         """score density
##         """
##         pass # XXX ???

    def _end_Hsp_qseq(self):
        """alignment string for the query
        """
        self._hsp.query = self._value

    def _end_Hsp_hseq(self):
        """alignment string for the database
        """
        self._hsp.sbjct = self._value

    def _end_Hsp_midline(self):
        """Formatting middle line as normally seen in BLAST report
        """
        self._hsp.match = self._value # do NOT strip spaces!
        assert len(self._hsp.match)==len(self._hsp.query)
        assert len(self._hsp.match)==len(self._hsp.sbjct)

    # Statistics
    def _end_Statistics_db_num(self):
        """number of sequences in the database
        """
        self._blast.num_sequences_in_database = int(self._value)

    def _end_Statistics_db_len(self):
        """number of letters in the database
        """
        self._blast.num_letters_in_database = int(self._value)

    def _end_Statistics_hsp_len(self):
        """the effective HSP length
        """
        self._blast.effective_hsp_length = int(self._value)

    def _end_Statistics_eff_space(self):
        """the effective search space
        """
        self._blast.effective_search_space = float(self._value)

    def _end_Statistics_kappa(self):
        """Karlin-Altschul parameter K
        """
        self._blast.ka_params = float(self._value)

    def _end_Statistics_lambda(self):
        """Karlin-Altschul parameter Lambda
        """
        self._blast.ka_params = (float(self._value),
                                 self._blast.ka_params)

    def _end_Statistics_entropy(self):
        """Karlin-Altschul parameter H
        """
        self._blast.ka_params = self._blast.ka_params + (float(self._value),)

def parse(handle, debug=0):
    """Returns an iterator a Blast record for each query.

    handle - file handle to and XML file to parse
    debug - integer, amount of debug information to print

    This is a generator function that returns multiple Blast records
    objects - one for each query sequence given to blast.  The file
    is read incrementally, returning complete records as they are read
    in.

    Should cope with new BLAST 2.2.14+ which gives a single XML file
    for mutliple query records.

    Should also cope with XML output from older versions BLAST which
    gave multiple XML files concatenated together (giving a single file
    which strictly speaking wasn't valid XML)."""
    from xml.parsers import expat
    BLOCK = 1024
    MARGIN = 10 # must be at least length of newline + XML start
    XML_START = "<?xml"

    text = handle.read(BLOCK)
    pending = ""

    if not text:
        #NO DATA FOUND!
        raise ValueError("Your XML file was empty")
    
    while text:
        #We are now starting a new XML file
        if not text.startswith(XML_START):
            raise ValueError("Your XML file did not start with %s... "
                             "but instead %s" \
                             % (XML_START, repr(text[:20])))

        expat_parser = expat.ParserCreate()
        blast_parser = BlastParser(debug)
        expat_parser.StartElementHandler = blast_parser.startElement
        expat_parser.EndElementHandler = blast_parser.endElement
        expat_parser.CharacterDataHandler = blast_parser.characters

        expat_parser.Parse(text, False)
        while blast_parser._records:
            record = blast_parser._records[0]
            blast_parser._records = blast_parser._records[1:]
            yield record

        while True:
            #Read in another block of the file...
            text, pending = pending + handle.read(BLOCK), ""
            if not text:
                #End of the file!
                expat_parser.Parse("", True) # End of XML record
                break

            #Now read a little bit more so we can check for the
            #start of another XML file...
            pending = handle.read(MARGIN)

            if (text+pending).find("\n" + XML_START) == -1:
                # Good - still dealing with the same XML file
                expat_parser.Parse(text, False)        
                while blast_parser._records:
                    yield blast_parser._records.pop(0)
            else:
                # This is output from pre 2.2.14 BLAST,
                # one XML file for each query!
                
                # Finish the old file:
                text, pending = (text+pending).split("\n" + XML_START,1)
                pending = XML_START + pending

                expat_parser.Parse(text, True) # End of XML record
                while blast_parser._records:
                    yield blast_parser._records.pop(0)
               
                #Now we are going to re-loop, reset the
                #parsers and start reading the next XML file
                text, pending = pending, ""
                break

        #this was added because it seems that the Jython expat parser
        #was adding records later then the Python one
        while blast_parser._records:
            yield blast_parser._records.pop(0)
            
        #At this point we have finished the first XML record.
        #If the file is from an old version of blast, it may
        #contain more XML records (check if text=="").
        assert pending==""
        assert len(blast_parser._records) == 0
        
    #We should have finished the file!
    assert text==""
    assert pending==""
    assert len(blast_parser._records) == 0
    
"""
Bio/Blast/Record.py
==============================================================================
"""

# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Record classes to hold BLAST output.

Classes:
Blast              Holds all the information from a blast search.
PSIBlast           Holds all the information from a psi-blast search.

Header             Holds information from the header.
Description        Holds information about one hit description.
Alignment          Holds information about one alignment hit.
HSP                Holds information about one HSP.
MultipleAlignment  Holds information about a multiple alignment.
DatabaseReport     Holds information from the database report.
Parameters         Holds information from the parameters.

"""
# XXX finish printable BLAST output

#from Bio.Align import Generic

class Header:
    """Saves information from a blast header.

    Members:
    application         The name of the BLAST flavor that generated this data.
    version             Version of blast used.
    date                Date this data was generated.
    reference           Reference for blast.

    query               Name of query sequence.
    query_letters       Number of letters in the query sequence.  (int)
    
    database            Name of the database.
    database_sequences  Number of sequences in the database.  (int)
    database_letters    Number of letters in the database.  (int)

    """
    def __init__(self):
        self.application = ''
        self.version = ''
        self.date = ''
        self.reference = ''

        self.query = ''
        self.query_letters = None

        self.database = ''
        self.database_sequences = None
        self.database_letters = None

class Description:
    """Stores information about one hit in the descriptions section.

    Members:
    title           Title of the hit.
    score           Number of bits.  (int)
    bits            Bit score. (float)
    e               E value.  (float)
    num_alignments  Number of alignments for the same subject.  (int)
    
    """
    def __init__(self):
        self.title = ''
        self.score = None
        self.bits = None
        self.e = None
        self.num_alignments = None
    def __str__(self):
        return "%-66s %5s  %s" % (self.title, self.score, self.e)

class Alignment:
    """Stores information about one hit in the alignments section.

    Members:
    title      Name.
    hit_id     Hit identifier. (str)
    hit_def    Hit definition. (str)
    length     Length.  (int)
    hsps       A list of HSP objects.

    """
    def __init__(self):
        self.title = ''
        self.hit_id = ''
        self.hit_def = ''
        self.length = None
        self.hsps = []
    def __str__(self):
        lines = self.title.split('\n')
        lines.append("Length = %s\n" % self.length)
        return '\n           '.join(lines)

class HSP:
    """Stores information about one hsp in an alignment hit.

    Members:
    score           BLAST score of hit.  (float)
    bits            Number of bits for that score.  (float)
    expect          Expect value.  (float)
    num_alignments  Number of alignments for same subject.  (int)
    identities      Number of identities (int) if using the XML parser.
                    Tuple of numer of identities/total aligned (int, int)
                    if using the (obsolete) plain text parser.
    positives       Number of positives (int) if using the XML parser.
                    Tuple of numer of positives/total aligned (int, int)
                    if using the (obsolete) plain text parser.
    gaps            Number of gaps (int) if using the XML parser.
                    Tuple of numer of gaps/total aligned (int, int) if
                    using the (obsolete) plain text parser.
    align_length    Length of the alignment. (int)
    strand          Tuple of (query, target) strand.
    frame           Tuple of 1 or 2 frame shifts, depending on the flavor.

    query           The query sequence.
    query_start     The start residue for the query sequence.  (1-based)
    query_end       The end residue for the query sequence.  (1-based)
    match           The match sequence.
    sbjct           The sbjct sequence.
    sbjct_start     The start residue for the sbjct sequence.  (1-based)
    sbjct_end       The end residue for the sbjct sequence.  (1-based)
    
    Not all flavors of BLAST return values for every attribute:
              score     expect     identities   positives    strand  frame
    BLASTP     X          X            X            X
    BLASTN     X          X            X            X          X
    BLASTX     X          X            X            X                  X
    TBLASTN    X          X            X            X                  X
    TBLASTX    X          X            X            X                 X/X

    Note: for BLASTX, the query sequence is shown as a protein sequence,
    but the numbering is based on the nucleotides.  Thus, the numbering
    is 3x larger than the number of amino acid residues.  A similar effect
    can be seen for the sbjct sequence in TBLASTN, and for both sequences
    in TBLASTX.

    Also, for negative frames, the sequence numbering starts from
    query_start and counts down.

    """
    def __init__(self):
        self.score = None
        self.bits = None
        self.expect = None
        self.num_alignments = None
        self.identities = (None, None)
        self.positives = (None, None)
        self.gaps = (None, None)
        self.align_length = None
        self.strand = (None, None)
        self.frame = ()
        
        self.query = ''
        self.query_start = None
        self.query_end = None
        self.match = ''
        self.sbjct = ''
        self.sbjct_start = None
        self.sbjct_end = None

    def __str__(self):
        lines = ["Score %i (%i bits), expectation %0.1e, alignment length %i" \
                 % (self.score, self.bits, self.expect, self.align_length)]
        if self.align_length < 50:
            lines.append("Query:%s %s %s" % (str(self.query_start).rjust(8),
                                       str(self.query),
                                       str(self.query_end)))
            lines.append("               %s" \
                         % (str(self.match)))
            lines.append("Sbjct:%s %s %s" % (str(self.sbjct_start).rjust(8),
                                       str(self.sbjct),
                                       str(self.sbjct_end)))
        else:
            lines.append("Query:%s %s...%s %s" \
                         % (str(self.query_start).rjust(8),
                            str(self.query)[:45],
                            str(self.query)[-3:],
                            str(self.query_end)))
            lines.append("               %s...%s" \
                         % (str(self.match)[:45],
                            str(self.match)[-3:]))
            lines.append("Sbjct:%s %s...%s %s" \
                         % (str(self.sbjct_start).rjust(8),
                            str(self.sbjct)[:45],
                            str(self.sbjct)[-3:],
                            str(self.sbjct_end)))
        return "\n".join(lines)

class MultipleAlignment:
    """Holds information about a multiple alignment.

    Members:
    alignment  A list of tuples (name, start residue, sequence, end residue).

    The start residue is 1-based.  It may be blank, if that sequence is
    not aligned in the multiple alignment.

    """
    def __init__(self):
        self.alignment = []

    def to_generic(self, alphabet):
        """Retrieve generic alignment object for the given alignment.

        Instead of the tuples, this returns an Alignment object from
        Bio.Align.Generic, through which you can manipulate and query
        the object.

        alphabet is the specified alphabet for the sequences in the code (for
        example IUPAC.IUPACProtein.

        Thanks to James Casbon for the code.
        """
        #TODO - Switch to new Bio.Align.MultipleSeqAlignment class?
        seq_parts = []
        seq_names = []
        parse_number = 0
        n = 0
        for name, start, seq, end in self.alignment:
            if name == 'QUERY': #QUERY is the first in each alignment block
                parse_number += 1
                n = 0

            if parse_number == 1: # create on first_parse, append on all others
                seq_parts.append(seq)
                seq_names.append(name)
            else:
                seq_parts[n] += seq
                n += 1

        generic = []#Generic.Alignment(alphabet)
        for (name,seq) in zip(seq_names,seq_parts):
            #generic.add_sequence(name, seq)
            generic.append((name, seq))

        return generic

class Round:
    """Holds information from a PSI-BLAST round.

    Members:
    number       Round number.  (int)
    reused_seqs  Sequences in model, found again.  List of Description objects.
    new_seqs     Sequences not found, or below threshold.  List of Description.
    alignments          A list of Alignment objects.
    multiple_alignment  A MultipleAlignment object.
    
    """
    def __init__(self):
        self.number = None
        self.reused_seqs = []
        self.new_seqs = []
        self.alignments = []
        self.multiple_alignment = None

class DatabaseReport:
    """Holds information about a database report.
    
    Members:
    database_name              List of database names.  (can have multiple dbs)
    num_letters_in_database    Number of letters in the database.  (int)
    num_sequences_in_database  List of number of sequences in the database.
    posted_date                List of the dates the databases were posted.
    ka_params                  A tuple of (lambda, k, h) values.  (floats)
    gapped                     # XXX this isn't set right!
    ka_params_gap              A tuple of (lambda, k, h) values.  (floats)

    """
    def __init__(self):
        self.database_name = []
        self.posted_date = []
        self.num_letters_in_database = []
        self.num_sequences_in_database = []
        self.ka_params = (None, None, None)
        self.gapped = 0
        self.ka_params_gap = (None, None, None)

class Parameters:
    """Holds information about the parameters.

    Members:
    matrix              Name of the matrix.
    gap_penalties       Tuple of (open, extend) penalties.  (floats)
    sc_match            Match score for nucleotide-nucleotide comparison
    sc_mismatch         Mismatch penalty for nucleotide-nucleotide comparison
    num_hits            Number of hits to the database.  (int)
    num_sequences       Number of sequences.  (int)
    num_good_extends    Number of extensions.  (int)
    num_seqs_better_e   Number of sequences better than e-value.  (int)
    hsps_no_gap         Number of HSP's better, without gapping.  (int)
    hsps_prelim_gapped  Number of HSP's gapped in prelim test.  (int)
    hsps_prelim_gapped_attemped  Number of HSP's attempted in prelim.  (int)
    hsps_gapped         Total number of HSP's gapped.  (int)
    query_length        Length of the query.  (int)
    query_id            Identifier of the query sequence. (str)
    database_length     Number of letters in the database.  (int)
    effective_hsp_length         Effective HSP length.  (int)
    effective_query_length       Effective length of query.  (int)
    effective_database_length    Effective length of database.  (int)
    effective_search_space       Effective search space.  (int)
    effective_search_space_used  Effective search space used.  (int)
    frameshift          Frameshift window.  Tuple of (int, float)
    threshold           Threshold.  (int)
    window_size         Window size.  (int)
    dropoff_1st_pass    Tuple of (score, bits).  (int, float)
    gap_x_dropoff       Tuple of (score, bits).  (int, float)
    gap_x_dropoff_final Tuple of (score, bits).  (int, float)
    gap_trigger         Tuple of (score, bits).  (int, float)
    blast_cutoff        Tuple of (score, bits).  (int, float)
    """
    def __init__(self):
        self.matrix = ''
        self.gap_penalties = (None, None)
        self.sc_match = None
        self.sc_mismatch = None
        self.num_hits = None
        self.num_sequences = None
        self.num_good_extends = None
        self.num_seqs_better_e = None
        self.hsps_no_gap = None
        self.hsps_prelim_gapped = None
        self.hsps_prelim_gapped_attemped = None
        self.hsps_gapped = None
        self.query_id = None
        self.query_length = None
        self.database_length = None
        self.effective_hsp_length = None
        self.effective_query_length = None
        self.effective_database_length = None
        self.effective_search_space = None
        self.effective_search_space_used = None
        self.frameshift = (None, None)
        self.threshold = None
        self.window_size = None
        self.dropoff_1st_pass = (None, None)
        self.gap_x_dropoff = (None, None)
        self.gap_x_dropoff_final = (None, None)
        self.gap_trigger = (None, None)
        self.blast_cutoff = (None, None)
    
class Blast(Header, DatabaseReport, Parameters):
    """Saves the results from a blast search.

    Members:
    descriptions        A list of Description objects.
    alignments          A list of Alignment objects.
    multiple_alignment  A MultipleAlignment object.
    + members inherited from base classes

    """
    def __init__(self):
        Header.__init__(self)
        DatabaseReport.__init__(self)
        Parameters.__init__(self)
        self.descriptions = []
        self.alignments = []
        self.multiple_alignment = None

class PSIBlast(Header, DatabaseReport, Parameters):
    """Saves the results from a blastpgp search.

    Members:
    rounds       A list of Round objects.
    converged    Whether the search converged.
    + members inherited from base classes

    """
    def __init__(self):
        Header.__init__(self)
        DatabaseReport.__init__(self)
        Parameters.__init__(self)
        self.rounds = []
        self.converged = 0

