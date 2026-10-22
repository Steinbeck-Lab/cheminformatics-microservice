import { useState } from "react";
import { ToolPageHead } from "@/components/common/Bento";
import { Link } from "react-router-dom";

const TermsOfService = () => {
  const [currentYear] = useState(new Date().getFullYear());

  return (
    <div className="pc-wrap pc-legal">
      {/* Content Area */}
      <div>
        <div>
          <ToolPageHead title="Terms of Service" description="Last updated: April 2025" />

          {/* Main Content */}
          <div className="glass pc-legal-card">
            <div className="pc-legal-body">
              <div className="prose dark:prose-invert max-w-none pc-prose">
                <h2 className="text-foreground ">1. Software License</h2>
                <p className="text-foreground ">
                  Cheminformatics Microservice is licensed under the MIT License - a short and
                  simple permissive license with conditions only requiring preservation of copyright
                  and license notices. Licensed works, modifications, and larger works may be
                  distributed under different terms and without source code.
                </p>

                <div className="bg-muted p-6 rounded-lg my-6 font-mono text-sm text-foreground ">
                  <p className="font-semibold mb-4">MIT License</p>
                  <p>Copyright (c) {currentYear} Cheminformatics Microservice</p>
                  <br />
                  <p>
                    Permission is hereby granted, free of charge, to any person obtaining a copy of
                    this software and associated documentation files (the "Software"), to deal in
                    the Software without restriction, including without limitation the rights to
                    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
                    the Software, and to permit persons to whom the Software is furnished to do so,
                    subject to the following conditions:
                  </p>
                  <br />
                  <p>
                    The above copyright notice and this permission notice shall be included in all
                    copies or substantial portions of the Software.
                  </p>
                  <br />
                  <p>
                    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
                    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
                    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
                    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
                    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
                    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
                  </p>
                </div>

                <h2 className="text-foreground ">2. Use of Service</h2>
                <p className="text-foreground ">
                  The Cheminformatics Microservice and related applications are provided to you for
                  research, educational, and non-commercial purposes. By using our service, you
                  agree to:
                </p>
                <ul className="text-foreground ">
                  <li>Use the service in compliance with all applicable laws and regulations</li>
                  <li>Not attempt to disrupt, degrade, or impair the service in any way</li>
                  <li>Not attempt to gain unauthorized access to any part of the service</li>
                  <li>Not use the service for any illegal or unauthorized purpose</li>
                </ul>

                <h2 className="text-foreground ">3. API Usage</h2>
                <p className="text-foreground ">
                  Our public API instance is provided on an "as is" and "as available" basis. While
                  we strive for high availability, we do not guarantee uninterrupted service. We
                  reserve the right to:
                </p>
                <ul className="text-foreground ">
                  <li>Impose rate limits on API requests to ensure fair usage</li>
                  <li>Modify, suspend, or discontinue the service with or without notice</li>
                  <li>Block access to users who violate these terms</li>
                </ul>

                <h2 className="text-foreground ">4. Third-Party Integrations</h2>
                <p className="text-foreground ">
                  The Cheminformatics Microservice integrates several third-party tools and
                  libraries, each with their own licenses and terms. Use of these components is
                  subject to their respective license terms. Key third-party components include:
                </p>
                <ul className="text-foreground ">
                  <li>RDKit (BSD license)</li>
                  <li>CDK (LGPL license)</li>
                  <li>OpenBabel (GPL license)</li>
                  <li>Ketcher (Apache 2.0 license)</li>
                </ul>

                <h2 className="text-foreground ">5. Data and Content</h2>
                <p className="text-foreground ">
                  You retain all rights to any data you submit to our service. We do not claim
                  ownership over your content. However, by using our service, you grant us a
                  non-exclusive license to use, process, and store your content solely for the
                  purpose of providing and improving our service.
                </p>

                <h2 className="text-foreground ">6. Disclaimer of Warranties</h2>
                <p className="text-foreground ">
                  THE SERVICE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
                  INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
                  PARTICULAR PURPOSE, AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
                  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES, OR OTHER LIABILITY, WHETHER IN
                  AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
                  WITH THE SERVICE OR THE USE OR OTHER DEALINGS IN THE SERVICE.
                </p>

                <h2 className="text-foreground ">7. Limitation of Liability</h2>
                <p className="text-foreground ">
                  TO THE MAXIMUM EXTENT PERMITTED BY LAW, IN NO EVENT WILL WE BE LIABLE FOR ANY
                  INDIRECT, INCIDENTAL, SPECIAL, CONSEQUENTIAL, OR PUNITIVE DAMAGES, INCLUDING
                  WITHOUT LIMITATION, LOSS OF PROFITS, DATA, USE, GOODWILL, OR OTHER INTANGIBLE
                  LOSSES, RESULTING FROM YOUR ACCESS TO OR USE OF OR INABILITY TO ACCESS OR USE THE
                  SERVICE.
                </p>

                <h2 className="text-foreground ">8. Changes to Terms</h2>
                <p className="text-foreground ">
                  We reserve the right to modify these terms at any time. We will provide notice of
                  significant changes by posting the new terms on the website. Your continued use of
                  the service after such modifications constitutes your acceptance of the revised
                  terms.
                </p>

                <h2 className="text-foreground ">9. Contact Us</h2>
                <p className="text-foreground ">
                  If you have any questions about these Terms, please contact the Steinbeck group at
                  the Friedrich Schiller University Jena. For more information, visit{" "}
                  <a
                    href="https://cheminf.uni-jena.de"
                    target="_blank"
                    rel="noopener noreferrer"
                    className="text-primary hover:underline"
                  >
                    https://cheminf.uni-jena.de
                  </a>
                  .
                </p>

                <div className="mt-8 pt-6 border-t border-border ">
                  <Link
                    to="/privacy"
                    className="text-primary hover:underline clay-interactive inline-block px-3 py-1 rounded-lg"
                  >
                    View our Privacy Policy
                  </Link>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default TermsOfService;
