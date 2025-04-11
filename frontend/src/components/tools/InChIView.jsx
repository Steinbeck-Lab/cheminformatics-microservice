import React, { useState, useEffect, useRef, useCallback } from "react";
import {
  HiOutlineClipboardCopy,
  HiOutlineCheck,
  HiOutlineExclamationCircle,
  HiOutlineInformationCircle,
  HiOutlineRefresh,
  HiOutlinePencil,
  HiOutlineX,
  HiOutlineCode,
  HiOutlineDocumentText,
  HiOutlineUpload,
  HiOutlineQuestionMarkCircle,
  HiOutlineAdjustments,
  HiOutlineLightningBolt,
  HiOutlineBeaker,
  HiOutlineTemplate,
} from "react-icons/hi";

// Import utility functions
import {
  INCHI_VERSIONS,
  loadInchiModule,
  convertMolfileToInchi,
  generateInchiKey,
  convertInchiToMolfile,
  convertAuxinfoToMolfile,
} from "../../utils/inchiUtils.js";

// Tooltip component for InChI options
const OptionTooltip = ({ content }) => (
  <div className="group relative inline-block">
    <HiOutlineQuestionMarkCircle className="h-4 w-4 ml-1 text-gray-500 dark:text-gray-400 inline-block align-text-bottom cursor-help" />
    <div className="absolute z-10 opacity-0 invisible group-hover:visible group-hover:opacity-100 transition-opacity duration-300 w-64 bg-white dark:bg-gray-800 p-2 rounded-md shadow-lg border border-gray-200 dark:border-gray-700 text-xs text-gray-700 dark:text-gray-300 -mt-1 left-6">
      {content}
    </div>
  </div>
);

// InChI options component
const InChIOptions = ({ onChange, inchiVersion, setInchiVersion }) => {
  const [mobileH, setMobileH] = useState(true);
  const [includeStereo, setIncludeStereo] = useState(true);
  const [stereoMode, setStereoMode] = useState("SAbs");
  const [recMet, setRecMet] = useState(false);
  const [ket, setKet] = useState(false);
  const [t15, setT15] = useState(false);
  const [polymers, setPolymers] = useState(false);
  const [NPZz, setNPZz] = useState(false);
  const [orgMet, setOrgMet] = useState(inchiVersion === "1.07.3-orgmet");
  const [showTautomerism, setShowTautomerism] = useState(true); // Add state for tautomerism visibility

  // Additional stereo options
  const [SUU, setSUU] = useState(false);
  const [SLUUD, setSLUUD] = useState(false);
  const [NEWPSOFF, setNEWPSOFF] = useState(false);

  // Additional tautomer options
  const [PT_06_00, setPT_06_00] = useState(false); // 1,3 heteroatom H-shift
  const [PT_13_00, setPT_13_00] = useState(false); // keten-inol exchange
  const [PT_16_00, setPT_16_00] = useState(false); // nitroso/oxime
  const [PT_18_00, setPT_18_00] = useState(false); // cyanic/iso-cyanic acids
  const [PT_22_00, setPT_22_00] = useState(false); // imine/imine
  const [PT_39_00, setPT_39_00] = useState(false); // nitrone/azoxy or Behrend rearrangement

  // Polymer detail options
  const [NoEdits, setNoEdits] = useState(false);
  const [FoldCRU, setFoldCRU] = useState(false);
  const [NoFrameShift, setNoFrameShift] = useState(false);

  // Generate and propagate options string when settings change
  useEffect(() => {
    const options = [];

    // Add options based on current settings
    if (!mobileH) options.push("FixedH");
    if (!includeStereo) options.push("SNon");

    if (includeStereo) {
      if (stereoMode === "SRel") options.push("SRel");
      if (stereoMode === "SRac") options.push("SRac");
      if (stereoMode === "SUCF") options.push("SUCF");

      // Additional stereo options
      if (SUU) options.push("SUU");
      if (SLUUD) options.push("SLUUD");
      if (NEWPSOFF) options.push("NEWPSOFF");
    }

    if (recMet) options.push("RecMet");
    if (ket) options.push("KET");
    if (t15) options.push("15T");

    // Additional tautomer options
    if (PT_06_00) options.push("PT_06_00");
    if (PT_13_00) options.push("PT_13_00");
    if (PT_16_00) options.push("PT_16_00");
    if (PT_18_00) options.push("PT_18_00");
    if (PT_22_00) options.push("PT_22_00");
    if (PT_39_00) options.push("PT_39_00");

    // Polymer options
    if (polymers) {
      options.push("Polymers");

      // Detailed polymer options
      if (NoEdits) options.push("NoEdits");
      if (FoldCRU) options.push("FoldCRU");
      if (NoFrameShift) options.push("NoFrameShift");
    }

    if (NPZz) options.push("NPZz");

    // orgMet is only available in the 1.07.3-orgmet version
    if (orgMet && inchiVersion === "1.07.3-orgmet") options.push("OrgMet");

    // Format the options string as required by the API
    const optionsString = options.map((opt) => `-${opt}`).join(" ");
    onChange(optionsString);
  }, [
    mobileH,
    includeStereo,
    stereoMode,
    recMet,
    ket,
    t15,
    polymers,
    NPZz,
    orgMet,
    SUU,
    SLUUD,
    NEWPSOFF,
    PT_06_00,
    PT_13_00,
    PT_16_00,
    PT_18_00,
    PT_22_00,
    PT_39_00,
    NoEdits,
    FoldCRU,
    NoFrameShift,
    inchiVersion,
    onChange,
  ]);

  // Update OrgMet option when InChI version changes
  useEffect(() => {
    if (inchiVersion === "1.07.3-orgmet") {
      setOrgMet(true); // Enable orgMet when switching to 1.07.3-orgmet
    } else if (orgMet) {
      setOrgMet(false); // Disable orgMet when switching away from 1.07.3-orgmet
    }
  }, [inchiVersion, orgMet]);

  // Update polymer-related options when polymers toggle changes
  useEffect(() => {
    if (!polymers) {
      setNoEdits(false);
      setFoldCRU(false);
      setNoFrameShift(false);
    }
  }, [polymers]);

  // Reset all options to default
  const resetOptions = () => {
    setMobileH(true);
    setIncludeStereo(true);
    setStereoMode("SAbs");
    setRecMet(false);
    setKet(false);
    setT15(false);
    setPolymers(false);
    setNPZz(false);
    setOrgMet(inchiVersion === "1.07.3-orgmet" ? true : false);
    setSUU(false);
    setSLUUD(false);
    setNEWPSOFF(false);
    setPT_06_00(false);
    setPT_13_00(false);
    setPT_16_00(false);
    setPT_18_00(false);
    setPT_22_00(false);
    setPT_39_00(false);
    setNoEdits(false);
    setFoldCRU(false);
    setNoFrameShift(false);
  };

  return (
    <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center">
          <div className="bg-purple-100 dark:bg-purple-900/50 p-2 rounded-lg mr-3">
            <HiOutlineCode className="h-5 w-5 text-purple-700 dark:text-purple-400" />
          </div>
          <h2 className="text-lg font-bold text-gray-800 dark:text-white">
            InChI Options
          </h2>
        </div>
        <button
          onClick={resetOptions}
          className="text-sm flex items-center text-gray-600 dark:text-gray-400 hover:text-purple-600 dark:hover:text-purple-400"
        >
          <HiOutlineRefresh className="mr-1 h-4 w-4" />
          Reset
        </button>
      </div>

      <div className="mb-4">
        <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
          InChI Version
        </label>
        <select
          className="w-full px-3 py-2 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500"
          value={inchiVersion}
          onChange={(e) => setInchiVersion(e.target.value)}
        >
          {Object.entries(INCHI_VERSIONS).map(([value, version]) => (
            <option key={value} value={value}>
              {version.label}
            </option>
          ))}
        </select>
      </div>

      <div className="space-y-4">
        <div className="flex items-center">
          <input
            id="mobileH"
            type="checkbox"
            checked={mobileH}
            onChange={(e) => setMobileH(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="mobileH"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Mobile H Perception
          </label>
          <OptionTooltip content="When checked, use mobile H perception. When unchecked, the FixedH option is used, which perceives fixed H positions." />
        </div>

        <div>
          <div className="flex items-center">
            <input
              id="includeStereo"
              type="checkbox"
              checked={includeStereo}
              onChange={(e) => setIncludeStereo(e.target.checked)}
              className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
            />
            <label
              htmlFor="includeStereo"
              className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
            >
              Include Stereo
            </label>
            <OptionTooltip content="Include stereochemical information in the InChI string. When unchecked, the SNon option is used, which ignores stereochemistry." />
          </div>

          {includeStereo && (
            <div className="ml-6 mt-2 grid grid-cols-2 gap-2">
              <div className="flex items-center">
                <input
                  id="stereo-abs"
                  type="radio"
                  name="stereoMode"
                  value="SAbs"
                  checked={stereoMode === "SAbs"}
                  onChange={() => setStereoMode("SAbs")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-abs"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Absolute
                </label>
                <OptionTooltip content="Use absolute stereochemistry (default option)" />
              </div>

              <div className="flex items-center">
                <input
                  id="stereo-rel"
                  type="radio"
                  name="stereoMode"
                  value="SRel"
                  checked={stereoMode === "SRel"}
                  onChange={() => setStereoMode("SRel")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-rel"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Relative
                </label>
                <OptionTooltip content="Use relative stereochemistry" />
              </div>

              <div className="flex items-center">
                <input
                  id="stereo-rac"
                  type="radio"
                  name="stereoMode"
                  value="SRac"
                  checked={stereoMode === "SRac"}
                  onChange={() => setStereoMode("SRac")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-rac"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Racemic
                </label>
                <OptionTooltip content="Use racemic stereochemistry" />
              </div>

              <div className="flex items-center">
                <input
                  id="stereo-ucf"
                  type="radio"
                  name="stereoMode"
                  value="SUCF"
                  checked={stereoMode === "SUCF"}
                  onChange={() => setStereoMode("SUCF")}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300"
                />
                <label
                  htmlFor="stereo-ucf"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  From chiral flag
                </label>
                <OptionTooltip content="Use stereochemistry as indicated by the chiral flag in the molfile" />
              </div>

              {/* Additional stereo options */}
              <div className="flex items-center col-span-2 mt-1">
                <input
                  id="SUU"
                  type="checkbox"
                  checked={SUU}
                  onChange={(e) => setSUU(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="SUU"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Always include omitted/undefined stereo (SUU)
                </label>
                <OptionTooltip content="Include all omitted/undefined stereo centers in the InChI string" />
              </div>

              <div className="flex items-center col-span-2">
                <input
                  id="SLUUD"
                  type="checkbox"
                  checked={SLUUD}
                  onChange={(e) => setSLUUD(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="SLUUD"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Different marks for unknown/undefined stereo (SLUUD)
                </label>
                <OptionTooltip content="Use different marking for unknown and undefined stereogenic centers" />
              </div>

              <div className="flex items-center col-span-2">
                <input
                  id="NEWPSOFF"
                  type="checkbox"
                  checked={NEWPSOFF}
                  onChange={(e) => setNEWPSOFF(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="NEWPSOFF"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Both ends of wedge point to stereocenters (NEWPSOFF)
                </label>
                <OptionTooltip content="Use old stereocenter recognition algorithm: a stereo bond should have the sharp end pointing to a stereocenter" />
              </div>
            </div>
          )}
        </div>

        <div className="flex items-center">
          <input
            id="recmet"
            type="checkbox"
            checked={recMet}
            onChange={(e) => setRecMet(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="recmet"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Include Bonds to Metal
          </label>
          <OptionTooltip content="Include metal-bonds in the InChI string, which are otherwise disconnected" />
        </div>

        {/* Tautomerism options section */}
        <div className="border-t border-gray-200 dark:border-gray-700 pt-3">
          <div className="flex items-center justify-between">
            <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2 flex items-center">
              Tautomerism Options
              <OptionTooltip content="Configure which types of tautomerism are considered when generating InChI" />
            </p>
            <button
              onClick={() => setShowTautomerism(!showTautomerism)}
              className="text-sm text-gray-600 dark:text-gray-400 hover:text-purple-600 dark:hover:text-purple-400"
            >
              {showTautomerism ? "Hide" : "Show"}
            </button>
          </div>

          {showTautomerism && (
            <div className="grid grid-cols-1 gap-2 ml-2">
              <div className="flex items-center">
                <input
                  id="ket"
                  type="checkbox"
                  checked={ket}
                  onChange={(e) => setKet(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="ket"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Keto-enol Tautomerism (KET)
                </label>
                <OptionTooltip content="Include keto-enol tautomerism when generating InChI" />
              </div>

              <div className="flex items-center">
                <input
                  id="15t"
                  type="checkbox"
                  checked={t15}
                  onChange={(e) => setT15(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="15t"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  1,5-tautomerism (15T)
                </label>
                <OptionTooltip content="Include 1,5-tautomerism when generating InChI" />
              </div>

              {/* Additional tautomer options available in 1.07 */}
              {(inchiVersion === "1.07.3" ||
                inchiVersion === "1.07.3-orgmet") && (
                <>
                  <div className="flex items-center">
                    <input
                      id="PT_06_00"
                      type="checkbox"
                      checked={PT_06_00}
                      onChange={(e) => setPT_06_00(e.target.checked)}
                      className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                    />
                    <label
                      htmlFor="PT_06_00"
                      className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                    >
                      1,3 heteroatom H-shift (PT_06_00)
                    </label>
                    <OptionTooltip content="Include 1,3 heteroatom H-shift tautomerism" />
                  </div>

                  <div className="flex items-center">
                    <input
                      id="PT_13_00"
                      type="checkbox"
                      checked={PT_13_00}
                      onChange={(e) => setPT_13_00(e.target.checked)}
                      className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                    />
                    <label
                      htmlFor="PT_13_00"
                      className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                    >
                      Keten-inol exchange (PT_13_00)
                    </label>
                    <OptionTooltip content="Include keten-inol exchange tautomerism" />
                  </div>

                  <div className="flex items-center">
                    <input
                      id="PT_16_00"
                      type="checkbox"
                      checked={PT_16_00}
                      onChange={(e) => setPT_16_00(e.target.checked)}
                      className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                    />
                    <label
                      htmlFor="PT_16_00"
                      className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                    >
                      Nitroso/oxime (PT_16_00)
                    </label>
                    <OptionTooltip content="Include nitroso/oxime tautomerism" />
                  </div>

                  <div className="flex items-center">
                    <input
                      id="PT_18_00"
                      type="checkbox"
                      checked={PT_18_00}
                      onChange={(e) => setPT_18_00(e.target.checked)}
                      className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                    />
                    <label
                      htmlFor="PT_18_00"
                      className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                    >
                      Cyanic/iso-cyanic acids (PT_18_00)
                    </label>
                    <OptionTooltip content="Include cyanic/iso-cyanic acids tautomerism" />
                  </div>

                  <div className="flex items-center">
                    <input
                      id="PT_22_00"
                      type="checkbox"
                      checked={PT_22_00}
                      onChange={(e) => setPT_22_00(e.target.checked)}
                      className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                    />
                    <label
                      htmlFor="PT_22_00"
                      className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                    >
                      Imine/imine (PT_22_00)
                    </label>
                    <OptionTooltip content="Include imine/imine tautomerism" />
                  </div>

                  <div className="flex items-center">
                    <input
                      id="PT_39_00"
                      type="checkbox"
                      checked={PT_39_00}
                      onChange={(e) => setPT_39_00(e.target.checked)}
                      className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                    />
                    <label
                      htmlFor="PT_39_00"
                      className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                    >
                      Nitrone/azoxy or Behrend rearrangement (PT_39_00)
                    </label>
                    <OptionTooltip content="Include nitrone/azoxy or Behrend rearrangement tautomerism" />
                  </div>
                </>
              )}
            </div>
          )}
        </div>

        {/* Polymer options section */}
        <div className="border-t border-gray-200 dark:border-gray-700 pt-3">
          <div className="flex items-center">
            <input
              id="polymers"
              type="checkbox"
              checked={polymers}
              onChange={(e) => setPolymers(e.target.checked)}
              className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
            />
            <label
              htmlFor="polymers"
              className="ml-2 block text-sm font-medium text-gray-700 dark:text-gray-300"
            >
              Treat Polymers
            </label>
            <OptionTooltip content="Process polymer structures according to InChI polymer extension" />
          </div>

          {polymers && (
            <div className="ml-6 mt-2 space-y-2">
              <div className="flex items-center">
                <input
                  id="NoEdits"
                  type="checkbox"
                  checked={NoEdits}
                  onChange={(e) => setNoEdits(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="NoEdits"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  No pre-edits of original polymer structure
                </label>
                <OptionTooltip content="Skip pre-editing of the original polymer structure" />
              </div>

              <div className="flex items-center">
                <input
                  id="FoldCRU"
                  type="checkbox"
                  checked={FoldCRU}
                  onChange={(e) => setFoldCRU(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="FoldCRU"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Enable CRU folding
                </label>
                <OptionTooltip content="Enable Constitutional Repeating Unit (CRU) folding for polymer representation" />
              </div>

              <div className="flex items-center">
                <input
                  id="NoFrameShift"
                  type="checkbox"
                  checked={NoFrameShift}
                  onChange={(e) => setNoFrameShift(e.target.checked)}
                  className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
                />
                <label
                  htmlFor="NoFrameShift"
                  className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
                >
                  Disable CRU frame shift
                </label>
                <OptionTooltip content="Disable frame shift for Constitutional Repeating Unit (CRU)" />
              </div>
            </div>
          )}
        </div>

        <div className="flex items-center">
          <input
            id="npzz"
            type="checkbox"
            checked={NPZz}
            onChange={(e) => setNPZz(e.target.checked)}
            className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
          />
          <label
            htmlFor="npzz"
            className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
          >
            Allow non-polymer Zz pseudoatoms
          </label>
          <OptionTooltip content="Allow Zz or * (star) pseudoatoms outside the polymer context" />
        </div>

        {inchiVersion === "1.07.3-orgmet" && (
          <div className="flex items-center border-t border-gray-200 dark:border-gray-700 pt-3">
            <input
              id="orgmet"
              type="checkbox"
              checked={orgMet}
              onChange={(e) => setOrgMet(e.target.checked)}
              className="h-4 w-4 text-green-600 focus:ring-green-500 border-gray-300 rounded"
            />
            <label
              htmlFor="orgmet"
              className="ml-2 block text-sm text-gray-700 dark:text-gray-300"
            >
              Molecular inorganics
            </label>
            <OptionTooltip content="Enable support for molecular inorganics (experimental feature in InChI 1.07.3-orgmet)" />
          </div>
        )}
      </div>
    </div>
  );
};

// Result block component for displaying InChI, InChIKey, and AuxInfo
const ResultBlock = ({ title, value, onCopy, copyState, icon }) => {
  return (
    <div className="mt-4 bg-white dark:bg-gray-800 p-4 rounded-lg border border-gray-200 dark:border-gray-700 shadow-lg hover:shadow-xl transition-shadow duration-300">
      <div className="flex justify-between items-center mb-2">
        <h3 className="text-sm font-medium text-gray-700 dark:text-gray-300 flex items-center">
          {icon && <span className="mr-2">{icon}</span>}
          {title}
        </h3>
        {value && (
          <button
            onClick={onCopy}
            className="inline-flex items-center px-2 py-1 text-xs text-green-700 dark:text-green-300 bg-green-100 dark:bg-green-900/30 rounded hover:bg-green-200 dark:hover:bg-green-800/50 transition-colors duration-200"
          >
            {copyState ? (
              <>
                <HiOutlineCheck className="w-3.5 h-3.5 mr-1" />
                Copied!
              </>
            ) : (
              <>
                <HiOutlineClipboardCopy className="w-3.5 h-3.5 mr-1" />
                Copy
              </>
            )}
          </button>
        )}
      </div>
      <div className="font-mono text-sm overflow-x-auto bg-gray-50 dark:bg-gray-900 p-3 rounded border border-gray-200 dark:border-gray-700 break-all max-h-24 overflow-y-auto">
        {value || (
          <span className="text-gray-400 dark:text-gray-500">
            No data generated yet
          </span>
        )}
      </div>
    </div>
  );
};

const InChIView = () => {
  // State variables
  const [inchi, setInchi] = useState("");
  const [auxInfo, setAuxInfo] = useState("");
  const [inchiKey, setInchiKey] = useState("");
  const [inputInchi, setInputInchi] = useState("");
  const [inputSmiles, setInputSmiles] = useState(""); // Added for SMILES input
  const [molfileContent, setMolfileContent] = useState("");
  const [logMessage, setLogMessage] = useState("");
  const [options, setOptions] = useState("");
  const [inchiVersion, setInchiVersion] = useState(
    Object.keys(INCHI_VERSIONS).find((key) => INCHI_VERSIONS[key].default) ||
      "1.07.3"
  );
  const [activeInputType, setActiveInputType] = useState("structure"); // 'structure', 'smiles', 'molfile', 'inchi'

  const [isLoading, setIsLoading] = useState(false);
  const [isEditorReady, setIsEditorReady] = useState(false);
  const [error, setError] = useState(null);
  const [copySuccess, setCopySuccess] = useState(false);
  const [showCopyModal, setShowCopyModal] = useState(false);
  const [copyModalText, setCopyModalText] = useState("");
  const [inchiModuleLoaded, setInchiModuleLoaded] = useState(false);
  const [retryAttempt, setRetryAttempt] = useState(0);
  const [autoGenerate, setAutoGenerate] = useState(false);

  // Refs
  const ketcherFrame = useRef(null);
  const copyTextRef = useRef(null);
  const fileInputRef = useRef(null);
  const messageHandlers = useRef({});
  const messageId = useRef(0);
  const previousOptionsRef = useRef("");

  // Reset error when needed
  useEffect(() => {
    if (error) setError(null);
  }, [inputInchi, inputSmiles, inchiVersion, activeInputType, error]);

  // Auto-select text in copy modal when it appears
  useEffect(() => {
    if (showCopyModal && copyTextRef.current) {
      copyTextRef.current.select();
    }
  }, [showCopyModal]);

  // Function to send messages to the iframe and wait for response
  const sendMessage = useCallback((type, payload = {}) => {
    return new Promise((resolve, reject) => {
      if (!ketcherFrame.current || !ketcherFrame.current.contentWindow) {
        reject(new Error("Ketcher frame not available"));
        return;
      }

      // Generate unique ID for this message
      const id = messageId.current++;

      // Store promise handlers
      messageHandlers.current[id] = { resolve, reject };

      // Send message to iframe
      const message = { id, type, payload };

      try {
        ketcherFrame.current.contentWindow.postMessage(message, "*");
      } catch (err) {
        delete messageHandlers.current[id];
        reject(new Error(`Failed to send message: ${err.message}`));
      }

      // Set timeout to reject if no response
      setTimeout(() => {
        if (messageHandlers.current[id]) {
          delete messageHandlers.current[id];
          reject(new Error("Timeout waiting for response"));
        }
      }, 10000); // 10 second timeout
    });
  }, []);

  // Function to communicate with ketcher through both direct access and messaging
  const executeKetcherCommand = useCallback(
    async (command, args = []) => {
      // First try direct access method
      try {
        if (
          ketcherFrame.current &&
          ketcherFrame.current.contentWindow.ketcher
        ) {
          const ketcher = ketcherFrame.current.contentWindow.ketcher;
          if (typeof ketcher[command] === "function") {
            return await ketcher[command](...args);
          }
        }
      } catch (directError) {
        console.debug(
          `Direct Ketcher access failed for ${command}`,
          directError
        );
        // Fall through to postMessage method
      }

      // Fall back to message passing
      try {
        return await sendMessage("ketcher-command", { command, args });
      } catch (msgError) {
        console.error(`Message passing failed for ${command}`, msgError);
        throw new Error(`Failed to execute ${command}: ${msgError.message}`);
      }
    },
    [sendMessage]
  );

  // Generate InChIKey from InChI
  const generateInChIKeyFromInChI = useCallback(
    async (inchiString) => {
      if (!inchiString) {
        console.warn("Cannot generate InChIKey: No InChI provided");
        return;
      }

      if (!inchiString.startsWith("InChI=")) {
        console.warn("Invalid InChI format for key generation:", inchiString);
        return;
      }

      try {
        console.log(
          "Generating InChIKey from InChI:",
          inchiString.substring(0, 40) + "..."
        );

        const result = await generateInchiKey(inchiString, inchiVersion);

        if (!result) {
          console.error("No result returned from InChIKey generation");
          return;
        }

        if (result.return_code === -1) {
          console.error("Failed to generate InChIKey:", result.message);
          return;
        }

        if (!result.inchikey) {
          console.error("InChIKey not found in result");
          return;
        }

        setInchiKey(result.inchikey);

        // If there's a message, add it to the log
        if (result.message) {
          setLogMessage((prevLog) =>
            prevLog ? `${prevLog}\n${result.message}` : result.message
          );
        }

        return result.inchikey;
      } catch (err) {
        console.error("Error generating InChIKey:", err);
        // Don't set an error message for this - it's a secondary operation
        return null;
      }
    },
    [inchiVersion]
  );

  // Generate InChI from molfile function
  const generateInChIFromMolfile = useCallback(
    async (molfile) => {
      try {
        if (!molfile || molfile.trim() === "") {
          throw new Error("No molecule data provided");
        }

        if (!inchiModuleLoaded) {
          throw new Error(
            "InChI module not loaded. Please wait or refresh the page."
          );
        }

        console.log(
          "Converting molfile to InChI with options:",
          options,
          "version:",
          inchiVersion
        );

        // Convert molfile to InChI
        const result = await convertMolfileToInchi(
          molfile,
          options,
          inchiVersion
        );

        if (!result) {
          throw new Error("No result returned from conversion function");
        }

        if (result.return_code === -1) {
          throw new Error(result.message || "Failed to generate InChI");
        }

        // Update state with results
        setInchi(result.inchi || "");
        setAuxInfo(result.auxinfo || "");

        // Append to log, not replace it
        setLogMessage(
          `InChI generated with options: ${options || "none"}\n${
            result.message || ""
          }`
        );

        // Generate InChIKey
        if (result.inchi) {
          await generateInChIKeyFromInChI(result.inchi);
        }

        return result;
      } catch (err) {
        console.error("Failed to generate InChI from molfile:", err);
        setError(`Failed to generate InChI: ${err.message}`);
        throw err;
      } finally {
        setIsLoading(false);
      }
    },
    [inchiModuleLoaded, options, inchiVersion, generateInChIKeyFromInChI]
  );

  // Generate InChI from Ketcher
  const generateInChI = useCallback(async () => {
    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!inchiModuleLoaded) {
      setError("InChI module not loaded. Please wait or refresh the page.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Generating InChI...");
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      // Get molfile from Ketcher
      const molfile = await executeKetcherCommand("getMolfile");

      if (!molfile || molfile.trim() === "") {
        setError("No structure drawn. Please draw a molecule first.");
        setIsLoading(false);
        return;
      }

      // Check if it's a reaction (not supported by InChI)
      if (molfile.includes("$RXN")) {
        setError("Reactions are not supported. Please draw a single molecule.");
        setIsLoading(false);
        return;
      }

      // Save the molfile for reference
      setMolfileContent(molfile);

      // Generate InChI from the molfile
      await generateInChIFromMolfile(molfile);
    } catch (err) {
      console.error("Failed to generate InChI:", err);
      setError(`Failed to generate InChI: ${err.message}`);
      setIsLoading(false);
    }
  }, [
    isEditorReady,
    inchiModuleLoaded,
    executeKetcherCommand,
    generateInChIFromMolfile,
  ]);

  // Set up message communication with the Ketcher iframe
  useEffect(() => {
    const handleMessage = (event) => {
      // Only accept messages from our iframe
      if (
        ketcherFrame.current &&
        event.source === ketcherFrame.current.contentWindow
      ) {
        const { id, type, status, data, error } = event.data;

        // Handle response to our message
        if (id && messageHandlers.current[id]) {
          const handler = messageHandlers.current[id];

          if (status === "success") {
            handler.resolve(data);
          } else if (status === "error") {
            handler.reject(new Error(error || "Unknown error"));
          }

          // Remove the handler after it's been used
          delete messageHandlers.current[id];
        }

        // Handle initialization message
        if (type === "ketcher-ready") {
          console.log("Ketcher ready message received");
          setIsEditorReady(true);
        }
      }
    };

    // Add listener for messages
    window.addEventListener("message", handleMessage);

    // Clean up on unmount
    return () => {
      window.removeEventListener("message", handleMessage);
      messageHandlers.current = {};
    };
  }, []);

  // Log option changes and auto-generate if enabled
  useEffect(() => {
    // Only proceed if options actually changed
    if (options !== previousOptionsRef.current) {
      if (previousOptionsRef.current && options) {
        // Log the option change
        const optionChangeLog = `InChI options changed: ${options || "none"}`;

        // Only append to existing log, don't replace it completely
        setLogMessage((prevLog) => {
          if (prevLog) {
            return `${optionChangeLog}\n${prevLog}`;
          }
          return optionChangeLog;
        });

        // Auto-generate if enabled and we have a valid structure
        if (autoGenerate && isEditorReady && molfileContent) {
          // Slight delay to allow UI to update first
          setTimeout(() => {
            generateInChI();
          }, 100);
        }
      }

      // Save current options for next comparison
      previousOptionsRef.current = options;
    }
  }, [options, autoGenerate, isEditorReady, molfileContent, generateInChI]);

  // Initialize Ketcher iframe communication
  const initializeKetcher = useCallback(() => {
    // Function to inject communication script into iframe
    const injectCommunicationScript = () => {
      try {
        const iframeWindow = ketcherFrame.current.contentWindow;
        const iframeDocument = iframeWindow.document;

        // Create a script element
        const script = iframeDocument.createElement("script");
        script.textContent = `
          // Set up message handler in the Ketcher iframe
          window.addEventListener('message', async function(event) {
            // Check source
            if (event.source !== window.parent) return;
            
            const { id, type, payload } = event.data;
            
            // Handle command execution
            if (type === 'ketcher-command') {
              try {
                const { command, args } = payload;
                
                if (!window.ketcher) {
                  throw new Error('Ketcher not initialized');
                }
                
                if (typeof window.ketcher[command] !== 'function') {
                  throw new Error(\`Command \${command} not available\`);
                }
                
                const result = await window.ketcher[command](...(args || []));
                event.source.postMessage({ id, status: 'success', data: result }, event.origin);
              } catch (error) {
                event.source.postMessage({ id, status: 'error', error: error.message }, event.origin);
              }
            }
          });
          
          // Wait for Ketcher to initialize and then notify parent
          const checkKetcher = () => {
            if (window.ketcher) {
              // Notify parent that Ketcher is ready
              window.parent.postMessage({ type: 'ketcher-ready' }, '*');
              console.log('Ketcher ready, notified parent');
            } else {
              // Check again in 100ms
              setTimeout(checkKetcher, 100);
            }
          };
          
          // Start checking for ketcher
          checkKetcher();
        `;

        // Add script to iframe
        iframeDocument.head.appendChild(script);
        console.log("Communication script injected into iframe");
      } catch (error) {
        console.error("Failed to inject communication script:", error);
      }
    };

    // Wait for iframe to load, then inject script
    if (ketcherFrame.current) {
      // Clear any previous onload
      ketcherFrame.current.onload = () => {
        console.log("Ketcher iframe loaded, injecting script");
        // Let the iframe load completely before injecting
        setTimeout(injectCommunicationScript, 500);
      };
    }
  }, []);

  // Use a more targeted approach for handling scrolling
  useEffect(() => {
    // Only focus on preventing initial load scroll, not restricting user scrolling
    let initialLoadComplete = false;

    const handleScroll = () => {
      // Only interfere with scrolling during the initial load
      if (!initialLoadComplete) {
        requestAnimationFrame(() => {
          // This is a gentle approach - it prevents jumps during loading
          // but doesn't continuously force the position
        });
      }
    };

    window.addEventListener("scroll", handleScroll, { passive: true });

    // After a short delay, allow normal scrolling
    const timeout = setTimeout(() => {
      initialLoadComplete = true;
    }, 1500); // 1.5 seconds should be enough for initial load

    return () => {
      window.removeEventListener("scroll", handleScroll);
      clearTimeout(timeout);
    };
  }, []);

  // Re-initialize when retry attempt changes without scrolling
  useEffect(() => {
    initializeKetcher();
  }, [retryAttempt, initializeKetcher]);

  // Enhanced copyToClipboard function with multiple fallback methods
  const copyToClipboard = async (text = null, type = "inchi") => {
    let textToCopy;

    switch (type) {
      case "inchi":
        textToCopy = text || inchi;
        break;
      case "key":
        textToCopy = text || inchiKey;
        break;
      case "auxinfo":
        textToCopy = text || auxInfo;
        break;
      default:
        textToCopy = text || inchi;
    }

    if (!textToCopy) {
      setError(
        `No ${type.toUpperCase()} to copy. Generate ${type.toUpperCase()} first.`
      );
      return;
    }

    // Try multiple clipboard copy methods in sequence
    try {
      // Method 1: Use the Clipboard API (modern browsers)
      if (navigator.clipboard && navigator.clipboard.writeText) {
        await navigator.clipboard.writeText(textToCopy);
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      }

      // Method 2: Use execCommand (older browsers)
      const textArea = document.createElement("textarea");
      textArea.value = textToCopy;

      // Make the textarea out of viewport
      textArea.style.position = "fixed";
      textArea.style.left = "-999999px";
      textArea.style.top = "-999999px";
      document.body.appendChild(textArea);

      // Select and copy
      textArea.focus();
      textArea.select();

      const successful = document.execCommand("copy");
      document.body.removeChild(textArea);

      if (successful) {
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
        return;
      } else {
        throw new Error("execCommand copy failed");
      }
    } catch (err) {
      console.error("Failed to copy text:", err);

      // Method 3: Show a modal with text to copy manually
      setCopyModalText(textToCopy);
      setShowCopyModal(true);
    }
  };

  // Check InChI module status
  useEffect(() => {
    const checkModuleStatus = async () => {
      try {
        setInchiModuleLoaded(false);

        // Try to load the module for the current version
        await loadInchiModule(inchiVersion)
          .then(() => {
            setInchiModuleLoaded(true);
            console.log(`InChI module ${inchiVersion} loaded successfully`);
          })
          .catch((err) => {
            console.error(`Failed to load InChI module ${inchiVersion}:`, err);
            setError(
              `Failed to load InChI module ${inchiVersion}. Please check your network connection.`
            );
          });
      } catch (err) {
        console.error("Error checking InChI module status:", err);
      }
    };

    checkModuleStatus();
  }, [inchiVersion]);

  // Load InChI into Ketcher
  const loadInChI = useCallback(async () => {
    if (!inputInchi.trim()) {
      setError("Please enter an InChI string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    if (!inputInchi.startsWith("InChI=")) {
      setError('Invalid InChI string. InChI should start with "InChI="');
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage(`Converting InChI to structure...`);

    try {
      // Call the convertInchiToMolfile function
      const result = await convertInchiToMolfile(inputInchi, "", inchiVersion);

      if (result.return_code === -1) {
        throw new Error(
          result.message || "Failed to convert InChI to structure"
        );
      }

      if (result.molfile) {
        // Load the molfile into Ketcher
        await executeKetcherCommand("setMolecule", [result.molfile]);
        console.log("InChI loaded successfully");

        // Store the molfile content
        setMolfileContent(result.molfile);

        // Update state
        setInchi(inputInchi);
        setLogMessage(
          `Structure loaded from InChI successfully.\n${result.message || ""}`
        );

        // Generate InChIKey
        generateInChIKeyFromInChI(inputInchi);
      } else {
        throw new Error("No structure data generated from InChI");
      }
    } catch (err) {
      console.error("Failed to load InChI:", err);
      setError(`Failed to convert InChI to structure: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  }, [
    inputInchi,
    isEditorReady,
    inchiVersion,
    executeKetcherCommand,
    generateInChIKeyFromInChI,
  ]);

  // Load SMILES into Ketcher
  const loadSmiles = useCallback(async () => {
    if (!inputSmiles.trim()) {
      setError("Please enter a SMILES string");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Loading SMILES structure...");

    // Clear previous results when loading new structure
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      // Load SMILES into Ketcher
      // Note: SMILES conversion is handled by Ketcher, not by InChI module
      await executeKetcherCommand("setMolecule", [inputSmiles]);

      // Get molfile from Ketcher
      const molfile = await executeKetcherCommand("getMolfile");

      // Store the molfile content
      setMolfileContent(molfile);

      setLogMessage(
        `SMILES loaded successfully. Use "Generate InChI" to create InChI.`
      );

      // Auto-generate InChI if enabled
      if (autoGenerate && molfile) {
        await generateInChIFromMolfile(molfile);
      }
    } catch (err) {
      console.error("Failed to load SMILES:", err);
      setError(`Failed to load structure from SMILES: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  }, [
    inputSmiles,
    isEditorReady,
    executeKetcherCommand,
    autoGenerate,
    generateInChIFromMolfile,
  ]);

  // Load molfile or AuxInfo into Ketcher
  const loadMolfile = useCallback(async () => {
    if (!molfileContent.trim()) {
      setError("Please enter or upload a molfile or AuxInfo");
      return;
    }

    if (!isEditorReady) {
      setError("Editor not ready. Please try again in a moment.");
      return;
    }

    setIsLoading(true);
    setError(null);
    setLogMessage("Processing input...");

    // Clear previous results when loading new structure
    setInchi("");
    setAuxInfo("");
    setInchiKey("");

    try {
      // Check if input is AuxInfo
      if (molfileContent.startsWith("AuxInfo=")) {
        // Handle AuxInfo
        setLogMessage("Converting AuxInfo to structure...");

        try {
          const result = await convertAuxinfoToMolfile(
            molfileContent,
            0,
            0,
            inchiVersion
          );

          if (result.molfile) {
            await executeKetcherCommand("setMolecule", [result.molfile]);
            setMolfileContent(result.molfile);
            setLogMessage(
              `AuxInfo converted to structure successfully.\n${
                result.message || ""
              }`
            );

            // Auto-generate InChI if enabled
            if (autoGenerate && result.molfile) {
              await generateInChIFromMolfile(result.molfile);
            }
          } else {
            throw new Error("Failed to convert AuxInfo to structure");
          }
        } catch (err) {
          setError(`AuxInfo conversion error: ${err.message}`);
        }
      }
      // Check if input is a molfile
      else if (
        molfileContent.includes("V2000") ||
        molfileContent.includes("V3000")
      ) {
        // Handle molfile
        await executeKetcherCommand("setMolecule", [molfileContent]);
        setLogMessage(
          "Molfile loaded successfully. Use 'Generate InChI' to create InChI."
        );

        // Auto-generate InChI if enabled
        if (autoGenerate) {
          await generateInChIFromMolfile(molfileContent);
        }
      } else {
        setError(
          "Invalid format. Please enter a valid Molfile or AuxInfo string"
        );
      }
    } catch (err) {
      console.error("Failed to load input:", err);
      setError(`Failed to load structure: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  }, [
    molfileContent,
    isEditorReady,
    inchiVersion,
    executeKetcherCommand,
    autoGenerate,
    generateInChIFromMolfile,
  ]);

  // Handle file upload for molfile
  const handleMolfileUpload = useCallback((event) => {
    const file = event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      setMolfileContent(e.target.result);
    };
    reader.readAsText(file);
  }, []);

  // Function to retry initialization
  const handleRetryInit = useCallback(() => {
    setIsEditorReady(false);
    setError(null);
    setRetryAttempt((prev) => prev + 1);
  }, []);

  // Clear the editor
  const clearEditor = useCallback(async () => {
    if (!isEditorReady) {
      console.warn("Editor not ready for clearing");
      return;
    }

    try {
      await executeKetcherCommand("setMolecule", [""]);
      setInchi("");
      setAuxInfo("");
      setInchiKey("");
      setMolfileContent("");
      setLogMessage("Editor cleared");
      console.log("Editor cleared successfully");
    } catch (err) {
      console.error("Failed to clear editor:", err);
      setError("Failed to clear the editor");
    }
  }, [isEditorReady, executeKetcherCommand]);

  // Handle auto-generate toggle
  const toggleAutoGenerate = useCallback(() => {
    setAutoGenerate((prev) => !prev);
    if (!autoGenerate) {
      setLogMessage((prev) =>
        prev
          ? `Auto-generate enabled. InChI will update automatically when options change.\n${prev}`
          : "Auto-generate enabled. InChI will update automatically when options change."
      );
    } else {
      setLogMessage((prev) =>
        prev
          ? `Auto-generate disabled. Use "Generate InChI" button to apply changes.\n${prev}`
          : 'Auto-generate disabled. Use "Generate InChI" button to apply changes.'
      );
    }
  }, [autoGenerate]);

  // Examples of common molecules with their InChI and SMILES
  const examples = [
    {
      name: "Ethanol",
      inchi: "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
      smiles: "CCO",
      description: "Alcohol",
    },
    {
      name: "Aspirin",
      inchi:
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
      smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
      description: "Pain reliever",
    },
    {
      name: "Caffeine",
      inchi:
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
      smiles: "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
      description: "Stimulant",
    },
    {
      name: "Ibuprofen",
      inchi:
        "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)",
      smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
      description: "Anti-inflammatory",
    },
  ];

  return (
    <div className="flex flex-col gap-6 p-4 md:p-6 bg-gradient-to-br from-green-50 to-emerald-100 dark:from-gray-900 dark:to-slate-900 min-h-screen">
      {/* Copy Modal - For fallback copying */}
      {showCopyModal && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4 animate-fadeIn">
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-2xl max-w-lg w-full">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-bold text-gray-900 dark:text-white">
                Copy Text
              </h3>
              <button
                onClick={() => setShowCopyModal(false)}
                className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
              >
                <HiOutlineX className="w-5 h-5" />
              </button>
            </div>
            <p className="mb-4 text-gray-700 dark:text-gray-300">
              Automatic copying failed. Please select and copy this text
              manually:
            </p>
            <div className="mb-4">
              <input
                type="text"
                ref={copyTextRef}
                value={copyModalText}
                readOnly
                className="w-full p-2 border border-gray-300 dark:border-gray-600 rounded font-mono text-sm bg-gray-50 dark:bg-gray-900"
                onClick={(e) => e.target.select()}
              />
            </div>
            <div className="flex justify-end gap-3">
              <button
                onClick={() => setShowCopyModal(false)}
                className="px-4 py-2 bg-gray-200 dark:bg-gray-700 text-gray-800 dark:text-gray-200 rounded hover:bg-gray-300 dark:hover:bg-gray-600"
              >
                Close
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Main content area */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
        {/* Left sidebar with controls */}
        <div className="lg:col-span-3 flex flex-col gap-4">
          {/* InChI Options */}
          <InChIOptions
            onChange={setOptions}
            inchiVersion={inchiVersion}
            setInchiVersion={setInchiVersion}
          />

          {/* Global Editor Controls */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
            <div className="flex items-center justify-between mb-4">
              <div className="flex items-center">
                <div className="bg-emerald-100 dark:bg-emerald-900/50 p-2 rounded-lg mr-3">
                  <HiOutlineAdjustments className="h-5 w-5 text-emerald-700 dark:text-emerald-400" />
                </div>
                <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                  Editor Controls
                </h2>
              </div>

              {/* Auto-generate toggle */}
              <div className="flex items-center">
                <button
                  onClick={toggleAutoGenerate}
                  className={`flex items-center text-sm font-medium ${
                    autoGenerate
                      ? "text-emerald-600 dark:text-emerald-400"
                      : "text-gray-500 dark:text-gray-400"
                  }`}
                  title="When enabled, InChI will be generated automatically when options change"
                >
                  <HiOutlineLightningBolt
                    className={`mr-1 h-4 w-4 ${
                      autoGenerate
                        ? "text-emerald-500 dark:text-emerald-400"
                        : "text-gray-400 dark:text-gray-500"
                    }`}
                  />
                  Auto
                </button>
              </div>
            </div>

            <div className="space-y-4">
              <button
                onClick={clearEditor}
                disabled={isLoading || !isEditorReady}
                className={`w-full px-4 py-2 rounded-lg text-sm font-medium transition-colors flex items-center justify-center ${
                  isLoading || !isEditorReady
                    ? "bg-gray-300 dark:bg-gray-700 text-gray-500 dark:text-gray-400 cursor-not-allowed"
                    : "bg-white dark:bg-gray-800 text-gray-700 dark:text-gray-300 hover:bg-gray-100 dark:hover:bg-gray-700 border border-gray-300 dark:border-gray-600"
                }`}
              >
                <HiOutlineRefresh className="h-4 w-4 mr-2" />
                Clear Editor
              </button>

              <button
                onClick={generateInChI}
                disabled={isLoading || !isEditorReady || !inchiModuleLoaded}
                className={`w-full px-4 py-2.5 rounded-lg text-sm font-medium transition-all duration-300 flex items-center justify-center ${
                  isLoading || !isEditorReady || !inchiModuleLoaded
                    ? "bg-gray-400 dark:bg-gray-600 text-white cursor-not-allowed"
                    : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800 text-white"
                }`}
              >
                <HiOutlineDocumentText className="h-4 w-4 mr-2" />
                {isLoading ? "Generating..." : "Generate InChI"}
              </button>
            </div>
          </div>

          {/* Tab-specific controls */}
          <div className="bg-white dark:bg-gray-800 p-6 rounded-xl shadow-lg border border-gray-100 dark:border-gray-700">
            <div className="flex items-center mb-4">
              <div className="bg-indigo-100 dark:bg-indigo-900/50 p-2 rounded-lg mr-3">
                <HiOutlinePencil className="h-5 w-5 text-indigo-700 dark:text-indigo-400" />
              </div>
              <h2 className="text-lg font-bold text-gray-800 dark:text-white">
                Structure Input
              </h2>
            </div>

            {/* Input Type Tabs */}
            <div className="mb-4">
              <div className="flex flex-wrap border-b border-gray-200 dark:border-gray-700">
                <button
                  className={`px-3 py-2 text-xs font-medium ${
                    activeInputType === "structure"
                      ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("structure")}
                >
                  <HiOutlineBeaker className="inline-block mr-1" />
                  Draw
                </button>
                <button
                  className={`px-3 py-2 text-xs font-medium ${
                    activeInputType === "smiles"
                      ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("smiles")}
                >
                  <HiOutlineTemplate className="inline-block mr-1" />
                  SMILES
                </button>
                <button
                  className={`px-3 py-2 text-xs font-medium ${
                    activeInputType === "molfile"
                      ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("molfile")}
                >
                  <HiOutlineDocumentText className="inline-block mr-1" />
                  Mol/Aux
                </button>
                <button
                  className={`px-3 py-2 text-xs font-medium ${
                    activeInputType === "inchi"
                      ? "text-green-600 dark:text-green-400 border-b-2 border-green-600 dark:border-green-400"
                      : "text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-white hover:border-b-2 hover:border-gray-300 dark:hover:border-gray-600"
                  }`}
                  onClick={() => setActiveInputType("inchi")}
                >
                  <HiOutlineCode className="inline-block mr-1" />
                  InChI
                </button>
              </div>
            </div>

            {/* SMILES Input Option */}
            {activeInputType === "smiles" && (
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="smiles-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter SMILES to Edit
                  </label>
                  <input
                    id="smiles-input"
                    type="text"
                    value={inputSmiles}
                    onChange={(e) => setInputSmiles(e.target.value)}
                    placeholder="e.g., CCO for ethanol"
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500 transition-all duration-200"
                  />
                </div>

                <button
                  onClick={loadSmiles}
                  disabled={isLoading || !isEditorReady || !inputSmiles.trim()}
                  className={`w-full relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-green-500 ${
                    isLoading || !isEditorReady || !inputSmiles.trim()
                      ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                      : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800"
                  }`}
                >
                  {isLoading
                    ? "Loading..."
                    : !isEditorReady
                    ? "Initializing..."
                    : "Load Structure"}
                </button>

                {/* Quick Examples */}
                <div className="mt-2">
                  <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                    Quick Examples:
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {examples.map((example, index) => (
                      <button
                        key={index}
                        onClick={() => setInputSmiles(example.smiles)}
                        className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-green-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                        title={`${example.name}: ${example.description}`}
                      >
                        {example.name}
                      </button>
                    ))}
                  </div>
                </div>
              </div>
            )}

            {/* Molfile/AuxInfo Input Option */}
            {activeInputType === "molfile" && (
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="molfile-input"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter Molfile or AuxInfo Content
                  </label>
                  <textarea
                    id="molfile-input"
                    value={molfileContent}
                    onChange={(e) => setMolfileContent(e.target.value)}
                    placeholder="Paste Molfile content here... or AuxInfo=1/0/N:1,2,3/E:(1,2)/rA:3nCCO/rB:s1;s2;/rC:;;;"
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500 transition-all duration-200"
                    rows={6}
                  />
                </div>

                <div className="flex flex-col space-y-3">
                  <div className="relative">
                    <input
                      type="file"
                      id="molfile-upload"
                      accept=".mol,.sdf,.mdl"
                      ref={fileInputRef}
                      onChange={handleMolfileUpload}
                      className="hidden"
                    />
                    <label
                      htmlFor="molfile-upload"
                      className="cursor-pointer flex items-center justify-center w-full px-4 py-2.5 bg-gray-100 dark:bg-gray-700 text-gray-800 dark:text-gray-200 rounded-lg hover:bg-gray-200 dark:hover:bg-gray-600 transition-colors duration-300"
                    >
                      <HiOutlineUpload className="w-5 h-5 mr-2" />
                      Upload File
                    </label>
                  </div>

                  <button
                    onClick={loadMolfile}
                    disabled={
                      isLoading || !isEditorReady || !molfileContent.trim()
                    }
                    className={`relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-green-500 ${
                      isLoading || !isEditorReady || !molfileContent.trim()
                        ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                        : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800"
                    }`}
                  >
                    <HiOutlineDocumentText className="mr-2 h-5 w-5" />
                    {isLoading ? "Loading..." : "Load Structure"}
                  </button>
                </div>
              </div>
            )}

            {/* InChI Input Option */}
            {activeInputType === "inchi" && (
              <div className="space-y-4">
                <div>
                  <label
                    htmlFor="inchi-input-structure"
                    className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-2"
                  >
                    Enter InChI to Convert
                  </label>
                  <textarea
                    id="inchi-input-structure"
                    value={inputInchi}
                    onChange={(e) => setInputInchi(e.target.value)}
                    placeholder="e.g., InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
                    className="w-full px-4 py-3 bg-white dark:bg-gray-700 border border-gray-300 dark:border-gray-600 rounded-lg text-gray-900 dark:text-white shadow-sm focus:ring-green-500 focus:border-green-500 transition-all duration-200"
                    rows={4}
                  />
                </div>

                <button
                  onClick={loadInChI}
                  disabled={
                    isLoading ||
                    !isEditorReady ||
                    !inputInchi.trim() ||
                    !inchiModuleLoaded
                  }
                  className={`w-full relative overflow-hidden px-4 py-2.5 rounded-lg text-white font-medium flex items-center justify-center transition-all duration-300 focus:outline-none focus:ring-2 focus:ring-offset-2 dark:focus:ring-offset-gray-800 focus:ring-green-500 ${
                    isLoading ||
                    !isEditorReady ||
                    !inputInchi.trim() ||
                    !inchiModuleLoaded
                      ? "bg-gray-400 dark:bg-gray-600 cursor-not-allowed"
                      : "bg-gradient-to-r from-green-600 to-emerald-700 hover:from-green-700 hover:to-emerald-800"
                  }`}
                >
                  {isLoading
                    ? "Loading..."
                    : !isEditorReady
                    ? "Initializing..."
                    : !inchiModuleLoaded
                    ? "Loading InChI Module..."
                    : "Convert to Structure"}
                </button>

                {/* Quick Examples */}
                <div className="mt-2">
                  <p className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-2">
                    Quick Examples:
                  </p>
                  <div className="flex flex-wrap gap-2">
                    {examples.map((example, index) => (
                      <button
                        key={index}
                        onClick={() => setInputInchi(example.inchi)}
                        className="inline-flex items-center px-3 py-1.5 border border-gray-300 dark:border-gray-600 shadow-sm text-xs font-medium rounded-full text-gray-700 dark:text-gray-300 bg-white dark:bg-gray-700 hover:bg-gray-50 dark:hover:bg-gray-600 focus:outline-none focus:ring-2 focus:ring-offset-2 focus:ring-green-500 dark:focus:ring-offset-gray-800 transition-all duration-200"
                        title={`${example.name}: ${example.description}`}
                      >
                        {example.name}
                      </button>
                    ))}
                  </div>
                </div>
              </div>
            )}

            {/* Drawing Instructions */}
            {activeInputType === "structure" && (
              <div className="text-sm text-gray-700 dark:text-gray-300">
                <p>
                  Use the Ketcher editor to draw your chemical structure, then
                  click "Generate InChI" to convert it.
                </p>
                <ul className="list-disc list-inside mt-2 space-y-1 text-gray-600 dark:text-gray-400">
                  <li>Click on atoms and bonds to create your structure</li>
                  <li>Use the toolbar for advanced drawing options</li>
                  <li>Clear the editor anytime to start over</li>
                </ul>
              </div>
            )}
          </div>

          {/* Status and Information */}
          <div
            className={`px-6 py-4 rounded-xl flex items-center justify-between ${
              isEditorReady && inchiModuleLoaded
                ? "bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800/40"
                : "bg-yellow-50 dark:bg-yellow-900/20 border border-yellow-200 dark:border-yellow-800/40"
            }`}
          >
            <div className="flex flex-col space-y-2">
              <div className="flex items-center">
                <div
                  className={`h-3 w-3 rounded-full mr-3 ${
                    isEditorReady
                      ? "bg-green-500 animate-pulse"
                      : "bg-yellow-500 animate-pulse"
                  }`}
                ></div>
                <span
                  className={`text-sm font-medium ${
                    isEditorReady
                      ? "text-green-800 dark:text-green-300"
                      : "text-yellow-800 dark:text-yellow-300"
                  }`}
                >
                  {isEditorReady ? "Editor Ready" : "Initializing Editor..."}
                </span>
              </div>

              <div className="flex items-center">
                <div
                  className={`h-3 w-3 rounded-full mr-3 ${
                    inchiModuleLoaded
                      ? "bg-green-500 animate-pulse"
                      : "bg-yellow-500 animate-pulse"
                  }`}
                ></div>
                <span
                  className={`text-sm font-medium ${
                    inchiModuleLoaded
                      ? "text-green-800 dark:text-green-300"
                      : "text-yellow-800 dark:text-yellow-300"
                  }`}
                >
                  {inchiModuleLoaded
                    ? "InChI Module Ready"
                    : "Loading InChI Module..."}
                </span>
              </div>
            </div>

            <div className="flex gap-3">
              {(!isEditorReady || !inchiModuleLoaded) && (
                <button
                  onClick={handleRetryInit}
                  className="inline-flex items-center px-3 py-1.5 rounded-lg text-sm font-medium transition-colors bg-yellow-100 dark:bg-yellow-900/30 text-yellow-700 dark:text-yellow-300 hover:bg-yellow-200 dark:hover:bg-yellow-900/50 border border-yellow-300 dark:border-yellow-700/50"
                >
                  <HiOutlineRefresh className="h-4 w-4 mr-1.5" />
                  Retry
                </button>
              )}
            </div>
          </div>

          {/* Information Box */}
          <div className="bg-gradient-to-br from-blue-50 to-indigo-50 dark:from-blue-900/30 dark:to-indigo-900/30 border border-blue-200 dark:border-blue-800/50 rounded-xl p-5 text-sm shadow-lg">
            <h4 className="font-bold text-blue-800 dark:text-blue-300 mb-3 flex items-center">
              <HiOutlineInformationCircle className="h-5 w-5 mr-2 text-blue-500 dark:text-blue-400" />
              About InChI
            </h4>
            <div className="space-y-3 text-gray-700 dark:text-gray-300">
              <p>
                The IUPAC International Chemical Identifier (InChI) is a textual
                identifier for chemical substances, designed to provide a
                standard way to encode molecular information.
              </p>
              <div>
                <h5 className="font-medium mb-1 text-gray-800 dark:text-gray-200">
                  Key Features:
                </h5>
                <ul className="list-disc list-inside space-y-1 pl-1 text-gray-600 dark:text-gray-400">
                  <li>
                    Canonical unique representation of chemical structures
                  </li>
                  <li>Hierarchical, layered information system</li>
                  <li>Machine-readable and human-readable format</li>
                  <li>Open-source, non-proprietary standard</li>
                  <li>
                    Facilitates data exchange between chemical information
                    systems
                  </li>
                </ul>
              </div>
            </div>
          </div>
        </div>

        {/* Main Editor Area */}
        <div className="lg:col-span-9 flex flex-col gap-4">
          {/* Error Display */}
          {error && (
            <div
              className="p-4 rounded-xl bg-red-50 dark:bg-red-900/30 text-red-700 dark:text-red-200 border border-red-300 dark:border-red-700/50 flex items-start shadow-lg animate-fadeIn"
              role="alert"
            >
              <HiOutlineExclamationCircle
                className="h-5 w-5 mr-3 flex-shrink-0 mt-0.5 text-red-500 dark:text-red-400"
                aria-hidden="true"
              />
              <p>{error}</p>
            </div>
          )}

          {/* Ketcher Editor Container */}
          <div
            className="bg-white dark:bg-gray-800 rounded-xl shadow-xl border border-gray-200 dark:border-gray-700 overflow-hidden transition-all duration-300 hover:shadow-2xl"
            style={{ height: "900px" }}
          >
            <iframe
              ref={ketcherFrame}
              src={`${process.env.PUBLIC_URL}/standalone/index.html`}
              title="Ketcher Editor"
              className="w-full h-full"
              frameBorder="0"
              onLoad={(e) => {
                // Don't prevent default - just let communication script be injected
                console.log("Ketcher iframe loaded");
              }}
            />
          </div>

          {/* Results */}
          <div className="bg-white dark:bg-gray-800 p-5 rounded-xl border border-gray-200 dark:border-gray-700 shadow-lg transition-all duration-300 hover:shadow-xl">
            <h3 className="text-lg font-bold text-gray-800 dark:text-white mb-3 flex items-center">
              <HiOutlineDocumentText className="h-5 w-5 mr-2 text-green-500 dark:text-green-400" />
              Results
            </h3>

            {/* InChI */}
            <ResultBlock
              title="InChI"
              value={inchi}
              onCopy={() => copyToClipboard(inchi, "inchi")}
              copyState={copySuccess}
              icon={
                <HiOutlineCode className="h-4 w-4 text-emerald-600 dark:text-emerald-500" />
              }
            />

            {/* InChIKey */}
            <ResultBlock
              title="InChIKey"
              value={inchiKey}
              onCopy={() => copyToClipboard(inchiKey, "key")}
              copyState={copySuccess}
              icon={
                <HiOutlineInformationCircle className="h-4 w-4 text-purple-600 dark:text-purple-500" />
              }
            />

            {/* AuxInfo */}
            <ResultBlock
              title="AuxInfo"
              value={auxInfo}
              onCopy={() => copyToClipboard(auxInfo, "auxinfo")}
              copyState={copySuccess}
              icon={
                <HiOutlineDocumentText className="h-4 w-4 text-blue-600 dark:text-blue-500" />
              }
            />

            {/* Log messages */}
            <div className="mt-4 p-3 bg-gray-50 dark:bg-gray-900 rounded-lg border border-gray-200 dark:border-gray-700 hover:border-green-300 dark:hover:border-green-700 transition-colors duration-300">
              <h4 className="text-sm font-medium text-gray-700 dark:text-gray-300 mb-1 flex items-center">
                <HiOutlineInformationCircle className="h-4 w-4 mr-1 text-gray-600 dark:text-gray-400" />
                Log
              </h4>
              <pre className="text-xs text-gray-600 dark:text-gray-400 font-mono whitespace-pre-wrap max-h-48 overflow-auto">
                {logMessage || "No log messages yet."}
              </pre>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default InChIView;
