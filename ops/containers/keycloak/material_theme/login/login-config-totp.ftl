<#import "template.ftl" as layout>
<@layout.registrationLayout displayInfo=true; section>
    <#if section = "title">
        ${msg("loginTotpTitle")?no_esc}
    <#elseif section = "header">
        ${msg("loginTotpTitle")?no_esc}
    <#elseif section = "form">
<ol id="kc-totp-settings">
    <li>
        <p>${msg("loginTotpStep1")}</p>
        </li>
    <li>
        <p>${msg("loginTotpStep2")}</p>
        <img src="data:image/png;base64, ${totp.totpSecretQrCode}" alt="Figure: Barcode"><br/>
        <span class="code">${totp.totpSecretEncoded}</span>
        </li>
    <li>
        <p>${msg("loginTotpStep3")}</p>
        </li>
    </ol>
    <form action="${url.loginAction}" class="form config-totp ${properties.kcFormClass!}" id="kc-totp-settings-form" method="post">
        <div class="${properties.kcFormGroupClass!} row">
            <div class="${properties.kcInputWrapperClass!} col-xs-12">
                <input type="text" id="totp" name="totp" autocomplete="off" class="form-control ${properties.kcInputClass!}" />
            </div>
            <input type="hidden" id="totpSecret" name="totpSecret" value="${totp.totpSecret}" />
        </div>
        <div class="row config-totp-button-container">
            <div class="col-xs-6 col-xs-offset-6 col-sm-4 col-sm-offset-8">
                <input class="btn btn-primary btn-flat btn-block ${properties.kcButtonClass!} ${properties.kcButtonPrimaryClass!} ${properties.kcButtonLargeClass!}" type="submit" value="${msg("doSubmit")}"/>
            </div>
        </div>
    </form>
    </#if>
</@layout.registrationLayout>
