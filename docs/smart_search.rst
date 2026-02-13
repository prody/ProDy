Smart Search
============

Type a ProDy function name.

- If it's an exact match (case-insensitive), you go directly to that function page.
- Otherwise, it falls back to normal search.

.. raw:: html

   <div style="max-width: 700px;">
     <form id="smartSearchForm">
       <input id="smartSearchInput"
              placeholder="Enter function name (e.g. parsePDB)"
              style="width:100%; padding:10px; font-size:16px;" />
       <button type="submit" style="margin-top:10px; padding:8px 12px;">Search</button>
     </form>
     <p id="smartSearchStatus" style="margin-top:10px; opacity:0.8;"></p>
   </div>
