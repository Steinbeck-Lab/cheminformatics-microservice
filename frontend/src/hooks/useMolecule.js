// Description: Custom hook for handling molecule data and operations
import { useState, useCallback, useEffect } from "react";
import { useAppContext } from "../context/AppContext";
import chemService from "../services/chemService";
import depictService from "../services/depictService";

/**
 * Custom hook for handling molecule data and operations
 * @param {string} initialSmiles - Initial SMILES string
 * @returns {Object} - Molecule state and functions
 */
const useMolecule = (initialSmiles = "") => {
  const [smiles, setSmiles] = useState(initialSmiles);
  const [moleculeData, setMoleculeData] = useState(null);
  const [isValid, setIsValid] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const { addRecentMolecule } = useAppContext();

  // Validate the SMILES string
  const validateMolecule = useCallback(
    async (smilesStr) => {
      if (!smilesStr) {
        setIsValid(false);
        setMoleculeData(null);
        return false;
      }

      setIsLoading(true);
      setError(null);

      try {
        // Use the errors endpoint to check if the SMILES is valid
        const result = await chemService.checkStructureErrors(smilesStr);

        // Check if the response contains error messages
        const isValidMolecule = result.messages.includes("No Errors Found");
        setIsValid(isValidMolecule);

        if (isValidMolecule) {
          // Add to recent molecules if it's valid
          addRecentMolecule({
            smiles: smilesStr,
            timestamp: new Date().toISOString(),
          });
        }

        return isValidMolecule;
      } catch (err) {
        setError(err.message);
        setIsValid(false);
        return false;
      } finally {
        setIsLoading(false);
      }
    },
    [addRecentMolecule]
  );

  // Update the SMILES string and validate
  const updateMolecule = useCallback(
    async (newSmiles) => {
      setSmiles(newSmiles);
      return validateMolecule(newSmiles);
    },
    [validateMolecule]
  );

  // Fetch additional data for the molecule
  const fetchMoleculeData = useCallback(async () => {
    if (!smiles || !isValid) return;

    setIsLoading(true);
    setError(null);

    try {
      // Fetch descriptors
      const descriptors = await chemService.calculateDescriptors(smiles);

      // Fetch NP-likeness score
      const npScore = await chemService.calculateNPLikeness(smiles);

      // Update molecule data
      setMoleculeData({
        smiles,
        descriptors,
        npScore,
        timestamp: new Date().toISOString(),
      });
    } catch (err) {
      setError(err.message);
    } finally {
      setIsLoading(false);
    }
  }, [smiles, isValid]);

  // Get standardized SMILES
  const getStandardizedSmiles = useCallback(async () => {
    if (!smiles) return null;

    setIsLoading(true);
    setError(null);

    try {
      const result = await chemService.checkStructureErrors(smiles, true);
      return result.standardized.smi;
    } catch (err) {
      setError(err.message);
      return null;
    } finally {
      setIsLoading(false);
    }
  }, [smiles]);

  // Get depiction URL
  const getDepictionUrl = useCallback(
    (options = {}) => {
      if (!smiles) return "";
      return depictService.get2DDepictionUrl(smiles, options);
    },
    [smiles]
  );

  // Get stereoisomers
  const getStereoisomers = useCallback(async () => {
    if (!smiles || !isValid) return [];

    setIsLoading(true);
    setError(null);

    try {
      return await chemService.generateStereoisomers(smiles);
    } catch (err) {
      setError(err.message);
      return [];
    } finally {
      setIsLoading(false);
    }
  }, [smiles, isValid]);

  // Clear current molecule data
  const clearMolecule = useCallback(() => {
    setSmiles("");
    setMoleculeData(null);
    setIsValid(false);
    setError(null);
  }, []);

  // Validate initial SMILES on mount
  useEffect(() => {
    if (initialSmiles) {
      validateMolecule(initialSmiles);
    }
  }, [initialSmiles, validateMolecule]);

  return {
    smiles,
    setSmiles: updateMolecule,
    moleculeData,
    isValid,
    isLoading,
    error,
    validateMolecule,
    fetchMoleculeData,
    getStandardizedSmiles,
    getDepictionUrl,
    getStereoisomers,
    clearMolecule,
  };
};

export default useMolecule;
