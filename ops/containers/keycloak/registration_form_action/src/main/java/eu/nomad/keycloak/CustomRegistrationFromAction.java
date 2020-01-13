package eu.nomad.keycloak;

import java.util.ArrayList;
import java.util.List;

import javax.ws.rs.core.MultivaluedMap;

import org.keycloak.Config.Scope;
import org.keycloak.authentication.FormAction;
import org.keycloak.authentication.FormActionFactory;
import org.keycloak.authentication.FormContext;
import org.keycloak.authentication.ValidationContext;
import org.keycloak.authentication.forms.RegistrationPage;
import org.keycloak.events.Errors;
import org.keycloak.forms.login.LoginFormsProvider;
import org.keycloak.models.AuthenticationExecutionModel.Requirement;
import org.keycloak.models.KeycloakSession;
import org.keycloak.models.KeycloakSessionFactory;
import org.keycloak.models.RealmModel;
import org.keycloak.models.UserModel;
import org.keycloak.models.utils.FormMessage;
import org.keycloak.provider.ProviderConfigProperty;

public class CustomRegistrationFromAction implements FormAction, FormActionFactory {

	private static final String PROVIDER_ID = "organization-field-validation-action";

	private static Requirement[] REQUIREMENT_CHOICES = { Requirement.REQUIRED, Requirement.DISABLED };

	@Override
	public void close() {
		// TODO Auto-generated method stub

	}

	@Override
	public void buildPage(FormContext arg0, LoginFormsProvider arg1) {
		// TODO Auto-generated method stub

	}

	@Override
	public boolean configuredFor(KeycloakSession arg0, RealmModel arg1, UserModel arg2) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean requiresUser() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void setRequiredActions(KeycloakSession arg0, RealmModel arg1, UserModel arg2) {
		// TODO Auto-generated method stub

	}

	@Override
	public void success(FormContext arg0) {
		// TODO Auto-generated method stub

    }

    private boolean validateUsername(String username) {
		return username.matches("^[a-zA-Z0-9_\\-\\.]+$");
    }

	@Override
	public void validate(ValidationContext context) {
		MultivaluedMap<String, String> formData = context.getHttpRequest().getDecodedFormParameters();
		List<FormMessage> errors = new ArrayList<>();

		String eventError = Errors.INVALID_REGISTRATION;

		if (!isBlank(formData.getFirst(RegistrationPage.FIELD_USERNAME))) {
			if (!this.validateUsername(formData.getFirst(RegistrationPage.FIELD_USERNAME))) {
				errors.add(new FormMessage(RegistrationPage.FIELD_USERNAME, "invalidUsernameMessage"));
			}
		}

		if (isBlank(formData.getFirst("user.attributes.affiliation"))) {
			errors.add(new FormMessage("user.attributes.affiliation", "missingAffiliationMessage"));
		}

		if (isBlank(formData.getFirst("user.attributes.affiliation_address"))) {
			errors.add(new FormMessage("user.attributes.affiliation_address", "missingAffiliationAddressMessage"));
		}

		if (errors.size() > 0) {
			context.error(eventError);
			context.validationError(formData, errors);
			return;

		} else {
			context.success();
		}
	}

	@Override
	public FormAction create(KeycloakSession arg0) {
		return this;
	}

	@Override
	public String getId() {
		return PROVIDER_ID;
	}

	@Override
	public void init(Scope arg0) {
	}

	@Override
	public void postInit(KeycloakSessionFactory arg0) {
	}

	@Override
	public String getDisplayType() {
		return "Signup Form Validation";
	}

	@Override
	public String getReferenceCategory() {
		return null;
	}

	@Override
	public Requirement[] getRequirementChoices() {
		return REQUIREMENT_CHOICES;
	}

	@Override
	public boolean isConfigurable() {
		return false;
	}

	@Override
	public boolean isUserSetupAllowed() {
		return false;
	}

	@Override
	public List<ProviderConfigProperty> getConfigProperties() {
		return null;
	}

	@Override
	public String getHelpText() {
		return "Validates Signup form fields.";
	}

	/**
	 * Check if string is blank (null or lenght is 0 or contains only white
	 * characters)
	 *
	 * @param s
	 *            to check
	 * @return true if string is blank
	 */
	public static boolean isBlank(String s) {
		return s == null || s.trim().length() == 0;
	}

}