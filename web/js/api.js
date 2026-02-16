// Data access layer: query Gaia DR3 and SIMBAD from the browser.
// Ported from src/data_access.py

const GAIA_TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync";
const SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync";

async function queryGaiaById(sourceId) {
  const query = `
    SELECT g.source_id, g.parallax, g.parallax_error,
           g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
           g.bp_rp, g.teff_gspphot,
           g.ag_gspphot, g.ebpminrp_gspphot,
           ap.mh_gspphot, ap.logg_gspphot, ap.lum_flame
    FROM gaiadr3.gaia_source AS g
    LEFT JOIN gaiadr3.astrophysical_parameters AS ap
      ON g.source_id = ap.source_id
    WHERE g.source_id = ${sourceId}
  `;

  const data = await tapQuery(GAIA_TAP_URL, query);
  if (!data || data.length === 0) return null;

  const row = data[0];
  return {
    source_id: String(row.source_id),
    parallax_mas: safeFloat(row.parallax),
    parallax_error_mas: safeFloat(row.parallax_error),
    phot_g_mean_mag: safeFloat(row.phot_g_mean_mag),
    phot_bp_mean_mag: safeFloat(row.phot_bp_mean_mag),
    phot_rp_mean_mag: safeFloat(row.phot_rp_mean_mag),
    bp_rp: safeFloat(row.bp_rp),
    ag_gspphot: safeFloat(row.ag_gspphot),
    ebpminrp_gspphot: safeFloat(row.ebpminrp_gspphot),
    feh: safeFloat(row.mh_gspphot),
    logg: safeFloat(row.logg_gspphot),
    teff_gspphot: safeFloat(row.teff_gspphot),
    lum_gspphot: safeFloat(row.lum_flame),
    is_cepheid: false,
    cepheid_period_days: null,
  };
}

async function resolveSimbadName(name) {
  // Step 1: query SIMBAD for the object's identifiers, looking for a Gaia DR3 ID
  const idQuery = `
    SELECT id FROM ident
    WHERE oidref = (SELECT oid FROM basic WHERE main_id = '${escapeSql(name)}'
                    UNION
                    SELECT oidref FROM ident WHERE id = '${escapeSql(name)}' LIMIT 1)
    AND id LIKE 'Gaia DR3%'
  `;

  try {
    const idData = await tapQuery(SIMBAD_TAP_URL, idQuery);
    if (idData && idData.length > 0) {
      const idStr = String(idData[0].id).trim();
      const match = idStr.match(/Gaia DR3\s+(\d+)/);
      if (match) return match[1];
    }
  } catch (e) {
    // Fall through to coordinate method
  }

  // Step 2: fallback -- get coordinates, do a Gaia cone search
  const coordQuery = `
    SELECT ra, dec FROM basic
    WHERE main_id = '${escapeSql(name)}'
    UNION
    SELECT b.ra, b.dec FROM basic b
    JOIN ident i ON b.oid = i.oidref
    WHERE i.id = '${escapeSql(name)}'
    LIMIT 1
  `;

  try {
    const coordData = await tapQuery(SIMBAD_TAP_URL, coordQuery);
    if (!coordData || coordData.length === 0) return null;

    const ra = coordData[0].ra;
    const dec = coordData[0].dec;

    const gaiaQuery = `
      SELECT TOP 1 source_id
      FROM gaiadr3.gaia_source
      WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', ${ra}, ${dec}, 0.005)) = 1
      ORDER BY phot_g_mean_mag ASC
    `;
    const gaiaData = await tapQuery(GAIA_TAP_URL, gaiaQuery);
    if (gaiaData && gaiaData.length > 0) {
      return String(gaiaData[0].source_id);
    }
  } catch (e) {
    // Resolution failed
  }

  return null;
}

async function tapQuery(url, adql) {
  const params = new URLSearchParams();
  params.append("REQUEST", "doQuery");
  params.append("LANG", "ADQL");
  params.append("FORMAT", "json");
  params.append("QUERY", adql);

  const response = await fetch(url, {
    method: "POST",
    body: params,
  });

  if (!response.ok) {
    throw new Error(`TAP query failed: ${response.status} ${response.statusText}`);
  }

  const json = await response.json();

  // TAP JSON format has metadata and data arrays
  if (json.data && json.metadata) {
    const columns = json.metadata.map(m => m.name);
    return json.data.map(row => {
      const obj = {};
      columns.forEach((col, i) => { obj[col] = row[i]; });
      return obj;
    });
  }

  // Some TAP services return VOTable-style JSON
  if (Array.isArray(json)) return json;

  return [];
}

function safeFloat(val) {
  if (val === null || val === undefined || val === "" || val === "NaN") return null;
  const n = Number(val);
  return isNaN(n) ? null : n;
}

function escapeSql(str) {
  return str.replace(/'/g, "''");
}
