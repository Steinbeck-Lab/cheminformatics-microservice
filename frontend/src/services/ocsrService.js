// OCSR Service using DECIMER
import api from "./api";

const OCSR_URL = "/ocsr";

export const processImage = async (image, reference = null, handDrawn = false) => {
  const formData = new FormData();
  formData.append("file", image);
  formData.append("hand_drawn", handDrawn);
  if (reference) {
    formData.append("reference", reference);
  }

  try {
    const response = await api.post(`${OCSR_URL}/process-upload`, formData, {
      headers: {
        "Content-Type": "multipart/form-data",
      },
    });
    return response.data;
  } catch (error) {
    throw new Error(`OCSR processing failed: ${error.message}`);
  }
};

export const processImageUrl = async (imageUrl, reference = null, handDrawn = false) => {
  try {
    const response = await api.post(`${OCSR_URL}/process`, {
      path: imageUrl,
      reference: reference || undefined,
      hand_drawn: handDrawn,
    });
    return response.data;
  } catch (error) {
    throw new Error(`OCSR processing failed: ${error.message}`);
  }
};

const ocsrService = {
  processImage,
  processImageUrl,
};

export default ocsrService;
