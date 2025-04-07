import React, { useState, useEffect } from "react";
import { motion } from "framer-motion";
import { Link } from "react-router-dom";

// Animation Variants
const pageVariants = {
  hidden: { opacity: 0 },
  visible: { opacity: 1, transition: { duration: 0.7, ease: "easeOut" } },
};

const headerVariants = {
  hidden: { opacity: 0, y: -20 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.7, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};

const contentVariants = {
  hidden: { opacity: 0, y: 30 },
  visible: {
    opacity: 1,
    y: 0,
    transition: { duration: 0.7, delay: 0.3, ease: [0.2, 0.65, 0.3, 0.9] },
  },
};

const TermsOfService = () => {
  const [currentYear] = useState(new Date().getFullYear());
  const [isLoaded, setIsLoaded] = useState(false);

  useEffect(() => {
    // Add a slight delay to ensure smooth animations after component mounts
    const timer = setTimeout(() => {
      setIsLoaded(true);
    }, 100);
    return () => clearTimeout(timer);
  }, []);

  return (
    <motion.div
      className="flex flex-col min-h-screen bg-slate-50 dark:bg-gray-900 text-slate-800 dark:text-slate-100 font-sans overflow-x-hidden"
      variants={pageVariants}
      initial="hidden"
      animate="visible"
    >
      {/* Background Effects */}
      <div className="absolute inset-0 -z-20 overflow-hidden pointer-events-none">
        <div className="absolute inset-0 dark:bg-gradient-to-br dark:from-slate-900 dark:via-indigo-950/50 dark:to-slate-950 opacity-0 dark:opacity-100 transition-opacity duration-500"></div>
        <div className="absolute inset-0 bg-gradient-to-br from-blue-50/80 via-white to-indigo-50/80 opacity-100 dark:opacity-0 transition-opacity duration-500"></div>
        <div className="absolute inset-0 w-full h-full bg-[radial-gradient(#e5e7eb_1px,transparent_1px)] dark:bg-[radial-gradient(#1e293b30_1px,transparent_1px)] [background-size:20px_20px] opacity-50 dark:opacity-30"></div>
      </div>

      {/* Content Area */}
      <div className="relative w-full mx-auto px-4 sm:px-6 lg:px-8 py-12 md:py-16 z-10 flex-grow">
        <div className="w-full max-w-4xl mx-auto">
          {/* Page Header */}
          <motion.div
            className="mb-12 text-center"
            variants={headerVariants}
            initial="hidden"
            animate={isLoaded ? "visible" : "hidden"}
          >
            <h1 className="text-4xl md:text-5xl font-bold mb-6">
              <span className="text-transparent bg-clip-text bg-gradient-to-r from-blue-600 to-indigo-600 dark:from-blue-400 dark:to-indigo-400">
                Terms of Service
              </span>
            </h1>
            <p className="text-lg text-slate-600 dark:text-slate-300">
              Last updated: April 2025
            </p>
          </motion.div>

          {/* Main Content */}
          <motion.div
            className="bg-white dark:bg-slate-800 rounded-xl shadow-lg border border-slate-200 dark:border-slate-700 overflow-hidden"
            variants={contentVariants}
            initial="hidden"
            animate={isLoaded ? "visible" : "hidden"}
          >
            <div className="p-8 md:p-10">
              <div className="prose prose-slate dark:prose-invert max-w-none">
                <h2 className="text-slate-800 dark:text-slate-200">
                  1. Software License
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  Cheminformatics Microservice is licensed under the MIT License
                  - a short and simple permissive license with conditions only
                  requiring preservation of copyright and license notices.
                  Licensed works, modifications, and larger works may be
                  distributed under different terms and without source code.
                </p>

                <div className="bg-slate-100 dark:bg-slate-700/50 p-6 rounded-lg my-6 font-mono text-sm text-slate-700 dark:text-slate-300">
                  <p className="font-semibold mb-4">MIT License</p>
                  <p>
                    Copyright (c) {currentYear} Cheminformatics Microservice
                  </p>
                  <br />
                  <p>
                    Permission is hereby granted, free of charge, to any person
                    obtaining a copy of this software and associated
                    documentation files (the "Software"), to deal in the
                    Software without restriction, including without limitation
                    the rights to use, copy, modify, merge, publish, distribute,
                    sublicense, and/or sell copies of the Software, and to
                    permit persons to whom the Software is furnished to do so,
                    subject to the following conditions:
                  </p>
                  <br />
                  <p>
                    The above copyright notice and this permission notice shall
                    be included in all copies or substantial portions of the
                    Software.
                  </p>
                  <br />
                  <p>
                    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
                    KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
                    WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
                    PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
                    OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
                    OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
                    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
                    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
                  </p>
                </div>

                <h2 className="text-slate-800 dark:text-slate-200">
                  2. Use of Service
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  The Cheminformatics Microservice and related applications are
                  provided to you for research, educational, and non-commercial
                  purposes. By using our service, you agree to:
                </p>
                <ul className="text-slate-700 dark:text-slate-300">
                  <li>
                    Use the service in compliance with all applicable laws and
                    regulations
                  </li>
                  <li>
                    Not attempt to disrupt, degrade, or impair the service in
                    any way
                  </li>
                  <li>
                    Not attempt to gain unauthorized access to any part of the
                    service
                  </li>
                  <li>
                    Not use the service for any illegal or unauthorized purpose
                  </li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">
                  3. API Usage
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  Our public API instance is provided on an "as is" and "as
                  available" basis. While we strive for high availability, we do
                  not guarantee uninterrupted service. We reserve the right to:
                </p>
                <ul className="text-slate-700 dark:text-slate-300">
                  <li>
                    Impose rate limits on API requests to ensure fair usage
                  </li>
                  <li>
                    Modify, suspend, or discontinue the service with or without
                    notice
                  </li>
                  <li>Block access to users who violate these terms</li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">
                  4. Third-Party Integrations
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  The Cheminformatics Microservice integrates several
                  third-party tools and libraries, each with their own licenses
                  and terms. Use of these components is subject to their
                  respective license terms. Key third-party components include:
                </p>
                <ul className="text-slate-700 dark:text-slate-300">
                  <li>RDKit (BSD license)</li>
                  <li>CDK (LGPL license)</li>
                  <li>OpenBabel (GPL license)</li>
                  <li>Ketcher (Apache 2.0 license)</li>
                </ul>

                <h2 className="text-slate-800 dark:text-slate-200">
                  5. Data and Content
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  You retain all rights to any data you submit to our service.
                  We do not claim ownership over your content. However, by using
                  our service, you grant us a non-exclusive license to use,
                  process, and store your content solely for the purpose of
                  providing and improving our service.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">
                  6. Disclaimer of Warranties
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  THE SERVICE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
                  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
                  WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
                  PURPOSE, AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
                  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER
                  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
                  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
                  SERVICE OR THE USE OR OTHER DEALINGS IN THE SERVICE.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">
                  7. Limitation of Liability
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  TO THE MAXIMUM EXTENT PERMITTED BY LAW, IN NO EVENT WILL WE BE
                  LIABLE FOR ANY INDIRECT, INCIDENTAL, SPECIAL, CONSEQUENTIAL,
                  OR PUNITIVE DAMAGES, INCLUDING WITHOUT LIMITATION, LOSS OF
                  PROFITS, DATA, USE, GOODWILL, OR OTHER INTANGIBLE LOSSES,
                  RESULTING FROM YOUR ACCESS TO OR USE OF OR INABILITY TO ACCESS
                  OR USE THE SERVICE.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">
                  8. Changes to Terms
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  We reserve the right to modify these terms at any time. We
                  will provide notice of significant changes by posting the new
                  terms on the website. Your continued use of the service after
                  such modifications constitutes your acceptance of the revised
                  terms.
                </p>

                <h2 className="text-slate-800 dark:text-slate-200">
                  9. Contact Us
                </h2>
                <p className="text-slate-700 dark:text-slate-300">
                  If you have any questions about these Terms, please contact
                  the Steinbeck group at the Friedrich Schiller University Jena.
                  For more information, visit{" "}
                  <a
                    href="https://cheminf.uni-jena.de"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-blue-600 dark:text-blue-400 hover:underline"
                  >
                    https://cheminf.uni-jena.de
                  </a>
                  .
                </p>

                <div className="mt-8 pt-6 border-t border-slate-200 dark:border-slate-700">
                  <Link
                    to="/privacy"
                    className="text-blue-600 dark:text-blue-400 hover:underline"
                  >
                    View our Privacy Policy
                  </Link>
                </div>
              </div>
            </div>
          </motion.div>
        </div>
      </div>
    </motion.div>
  );
};

export default TermsOfService;
